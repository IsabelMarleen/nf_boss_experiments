library(ggplot2)
library(argparse)
library(dplyr)
library(patchwork)
library(stringr)
library(tidyr)

theme_set(theme_minimal())




get_arguments <- function() {
  parser <- argparse::ArgumentParser()
  parser$add_argument('--input_cov_A', required = TRUE)
  parser$add_argument('--input_unb_A', required = TRUE)
parser$add_argument('--input_cov_B', required = TRUE)
  parser$add_argument('--input_unb_B', required = TRUE)
  parser$add_argument('--dump_time', required = TRUE)
  parser$add_argument('--pores', required = TRUE)
  parser$add_argument('--analysed_log', required = TRUE)
  parser$add_argument('--analysed_log2', required = TRUE)
  parser$add_argument('--output', required = TRUE)
  args<- parser$parse_args(commandArgs(trailingOnly = TRUE))
  return(args)
}

visualise_simulations <- function(ptol, nrow) {
  # get the arguments for input and output files
  args <- get_arguments()
  input_cov_a <- args$input_cov_A
  input_unb_a <- args$input_unb_A
  input_cov_b <- args$input_cov_B
  input_unb_b <- args$input_unb_B
  dump_time <- as.numeric(args$dump_time)
  pores <- as.numeric(args$pores)
  genome_size <- as.numeric(args$genome_size)
  log_a <- args$analysed_log
  log_b <- args$analysed_log2
  output <- args$output

  dump_time1 <- dump_time
  seq_speed  <- 400
  name_a  <- "HG003"
  name_b <- "HG004"

  # GET COVERAGE ----
  get_dump_time <- function(log){
    # load the analysed log
    log_file <- read.csv(log)

    # extract times for dumps
    log_time  <- log_file %>%
      group_by(cond, dump) %>%
      summarise(time=last(time)) %>%
      rename(seq_time = time, time = dump)
    
    return(log_time)
  }

  load_cov <- function(log_time, input_cov){
    # load the coverage data
    cov <- read.csv(input_cov) %>%
      left_join(log_time) %>% 
      mutate(seq_time = if_else(is.na(seq_time), 0, seq_time))

    return(cov)
  }

  # load the coverage data a
  log_time_a <- get_dump_time(log_a)
  cov_a <- load_cov(log_time_a, input_cov_a)
  
  # load the coverage data b
  log_time_b <- get_dump_time(log_b)
  cov_b <- load_cov(log_time_b, input_cov_b)

  # Combine coverage data frames
  cov <- bind_rows(name_a=cov_a, name_b=cov_b, .id="sim") %>%
    mutate(seq_time = seq_time/pores/60/seq_speed,
            evenness = if_else(seq_time==0, NA, evenness),
            sim = if_else(sim=="name_a", name_a, name_b))  # convert time units to minutes
  
  # get the ordering of the OTUs
  otu_order <- cov %>%
    filter(cond == 'control', time == max(time)) %>%
    arrange(desc(mean_coverage)) %>%
    distinct(otu) %>%
    pull(otu)

  # reorder factor levels
  cov$otu <- factor(cov$otu, levels=otu_order)
  
  cov <- cov %>% 
        group_by(cond,bc,otu)

  if ("bc" %in% colnames(cov)){
  bc_names  <- cov %>%
    filter(bc != "0.0")%>%
    arrange(bc)%>%
    pull(bc)%>%
    unique()
  
  # reorder factor levels
  cov$bc <- factor(cov$bc, levels=bc_names)

  } else {
    cov <- cov %>% 
        group_by(cond,otu)
  }


  # GET UNBLOCKS ----  
  # load the unblocking data and reorder the factor levels
  unb <- bind_rows(name_a=read.csv(input_unb_a), name_b=read.csv(input_unb_b), .id="sim")%>%
    left_join(bind_rows(name_a=log_time_a, name_b=log_time_b, .id="sim")) %>%
    mutate(seq_time = if_else(is.na(seq_time), 0, seq_time),
            sim = if_else(sim=="name_a", name_a, name_b))%>%
    mutate(orig_time = seq_time) %>%
    mutate(seq_time = if_else(is.na(seq_time), 0, seq_time/pores/60/seq_speed))%>%
    arrange(cond, time)
  unb$otu <- factor(unb$otu, levels=otu_order)

    if ("bc" %in% colnames(cov)){
    unb <- unb %>% 
      group_by(cond, otu, sim)%>%
      filter(bc=="0")%>%
      uncount(length(bc_names))%>%
      mutate(bc = bc_names)%>%
      bind_rows(filter(unb, bc!="0"))%>%
      mutate(orig_time = if_else(is.na(orig_time), 0, orig_time),
       unb_ratio = if_else(is.na(unb_ratio), 0, unb_ratio)) %>%
       arrange(cond,time,bc)%>%print()

    # reorder factor levels
    unb$bc <- factor(unb$bc, levels=bc_names)

    unb <- unb %>% 
      group_by(cond, bc, otu)
  }
  else{
    unb <- unb %>%
        group_by(cond,otu)
  }

  unb <- unb %>%
    rename("total_noncum" = total, "base_total_noncum" = base_total, "unb_noncum" = unb) %>%
    mutate(total = cumsum(total_noncum), 
      base_total = cumsum(as.numeric(base_total_noncum)),
      unb = cumsum(unb_noncum),
      seq_time_noncum = seq_time - lag(seq_time, default = first(seq_time)),
      unb_ratio = if_else(!is.na(unb/total), unb/total, 0),
      unb_ratio_noncum = if_else(!is.na(unb_noncum/total_noncum), unb_noncum/total_noncum, 0),
      total_rate = total/seq_time_noncum)

  unb_cond <- unb %>% 
    filter(cond == "boss")


  # create plots with new attempt at scaling time units to real sequencing time
   if ("bc" %in% colnames(cov)){
     col <- "bc"
  } else {
     col <- "otu"
  }

  unb_plot_seq <- ggplot(
      data=unb_cond,
      mapping=aes(x=seq_time, y=unb_ratio)) +
    geom_line(aes(color=.data[[col]])) +
    geom_point(aes(shape=.data[[col]],fill=.data[[col]], size=.data[[col]]), color="white", stroke=0.005, alpha=0.8) +
    scale_color_manual(values=ptol, guide = "none") +
    scale_fill_manual(values=ptol, guide = "none") +
    scale_size_manual(values=c(1.8,1.3,1.1), guide = "none")+
    scale_shape_manual(values=c(21,22,23), guide = "none")+
    ylab("cum. rejection rate") +
    xlim(0, max(unb_cond$seq_time))+
    facet_wrap(~otu, scales="free_y", nrow = nrow) +
    ylim(0,1)+
    theme(legend.position = "none")
  # unb_plot

  unb_plot_frag <- ggplot(
      data=unb_cond,
      # mapping=aes(x=seq_time, y=unb_ratio, color=bc)) +
      mapping=aes(x=total, y=unb_ratio)) +
    geom_line(aes(color=.data[[col]])) +
    geom_point(aes(shape=.data[[col]],fill=.data[[col]], size=.data[[col]]), color="white", stroke=0.005, alpha=0.8) +
    scale_color_manual(values=ptol, guide = "none") +
    scale_fill_manual(values=ptol, guide = "none") +
    # scale_linewidth_manual(values=c(1.4,1,0.6))+
    scale_size_manual(values=c(1.8,1.3,1.1), guide = "none")+
    scale_shape_manual(values=c(21,22,23), guide = "none")+
    ylab("cum. rejection rate") +
    xlim(0, max(unb_cond$total))+
    facet_wrap(~otu, scales="free_y", nrow = nrow) +
    ylim(0,1)+
    theme(legend.position = "none")

  unb_noncum_plot_seq <- ggplot(
      data=unb_cond,
      mapping=aes(x=seq_time, y=unb_ratio_noncum, color=.data[[col]])) +
    geom_line() + #linewidth=1,alpha=0.6
    geom_point(aes(shape=.data[[col]],fill=.data[[col]], size=.data[[col]]), color="white", stroke=0.005, alpha=0.8) +
    xlim(0, max(unb_cond$seq_time))+
    scale_color_manual(values=ptol) +
    scale_fill_manual(values=ptol) +
    scale_size_manual(values=c(2,1.4,1))+
    scale_shape_manual(values=c(21,22,23))+
    facet_wrap(~otu, scales="free_y", nrow = nrow) +
    ylab("rejection rate") +
    theme(legend.position = "none")
  # unb_plot

  nreads_seq <- ggplot(
      data=unb,
            mapping=aes(x=seq_time, y=total, linetype=cond, linewidth=.data[[col]], colour=.data[[col]])) +
    geom_line() +
    facet_wrap(~otu, scales="free_y", nrow = nrow) +
    xlim(0, max(unb_cond$seq_time))+
    scale_color_manual(values=ptol)+#, guide = "none") +
    scale_linewidth_manual(values=c(1,0.8,0.6))+
    ylab("# fragments")
  # nreads

  # nreads_noncum_seq <- ggplot(
  #     data=mutate(unb_noncum, total = if_else(time == 0, NA, total)),
  #     mapping=aes(x=seq_time, y=total, linetype=cond, colour=sim)) +
  #   geom_line(linewidth=1) +
  #   facet_wrap(~otu, scales="free_y", nrow = nrow) +
  #   scale_color_manual(values=ptol, guide = "none") +
  #   ylab("# reads")
  # # nreads

  meanc_seq <- ggplot(
      data=cov,
            mapping=aes(x=seq_time, y=mean_coverage, linetype=cond, colour=.data[[col]], linewidth=.data[[col]], alpha=.data[[col]])) +
    geom_line() +
    facet_wrap(~otu, scales="free_y", nrow = nrow) +
    xlim(0, max(unb_cond$seq_time))+
    ylim(0, max(pull(slice_min(filter(cov, seq_time > max(unb_cond$seq_time)), order_by=time), mean_coverage)))+
    scale_color_manual(values=ptol)+ #, guide = "none") +
    scale_linewidth_manual(values=c(1,0.8,0.6))+
    scale_alpha_manual(values=c(1,0.9,0.8), guide = "none")+
    ylab("mean coverage")+
    theme(legend.position = "none")
  # meanc

  lowc_seq <- ggplot(
      data=cov,
            mapping=aes(x=seq_time, y=low_coverage_prop, linetype=cond, colour=.data[[col]], linewidth=.data[[col]])) +
    geom_line() +
    facet_wrap(~otu, scales="free_y", nrow = nrow) +
    scale_color_manual(values=ptol) +
    xlim(0, max(unb_cond$seq_time))+
    scale_linewidth_manual(values=c(1,0.8,0.6))+
    ylab("prop. sites at <5x")+
    theme(legend.position = "none")
  # lowc

  evn_seq <- ggplot(
      data=cov,
            mapping=aes(x=seq_time, y=evenness, linetype=cond, colour=.data[[col]], linewidth=.data[[col]], alpha=cond)) +
    geom_line() +
    facet_wrap(~otu, scales="free_y", nrow = nrow) +
    scale_color_manual(values=ptol) +
    xlim(0, max(unb_cond$seq_time))+
    scale_linewidth_manual(values=c(1,0.8,0.6))+
    scale_alpha_manual(values=c(0.6,1), guide = "none")+
    ylab("evenness")+
    theme(legend.position = "none")
  # evenness of coverage

  layout_seqt <- wrap_plots(list(unb_plot_seq, unb_plot_frag, nreads_seq,lowc_seq, meanc_seq, evn_seq), axes='keep', axis_titles="keep")+
    plot_annotation(tag_levels = "A") +
    plot_layout(
      nrow = 3,
      ncol = 2,
      guides="collect"
      ) &
    xlab("seq. time (pseudominutes)") &
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      # strip.text.x = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    )

  layout_seqt[[2]] <- layout_seqt[[2]] + xlab("# fragments")
  ggsave(output, layout_seqt, w=8, h=7)

}


ptol <- c("#332288","#DDCC77","#117733","#88CCEE","#882255","#44AA99","#999933","#AA4499")#"#CC6677",
nrow <- 1
visualise_simulations(ptol, nrow)
