#!/usr/bin/env Rscript

library(ggplot2)
library(argparse)
library(dplyr)
library(patchwork)
library(stringr)
library(tidyr)

theme_set(theme_minimal())


get_arguments <- function() {
  parser <- argparse::ArgumentParser()
  parser$add_argument('--input_cov', required = TRUE)
  parser$add_argument('--input_unb', required = TRUE)
  parser$add_argument('--dump_time', required = TRUE)
  parser$add_argument('--pores', required = TRUE)
  parser$add_argument('--analysed_log', required = TRUE)
  parser$add_argument('--output', required = TRUE)
  args<- parser$parse_args(commandArgs(trailingOnly = TRUE))
  return(args)
}


visualise_simulation <- function(ptol, nrow) {
  # ARGUMENTS ----
  # get the arguments for input and output files
  args <- get_arguments()
  input_cov <- args$input_cov
  input_unb <- args$input_unb
  dump_time <- as.numeric(args$dump_time)
  pores <- as.numeric(args$pores)
  log_path <- args$analysed_log
  output <- args$output
  seq_speed  <- 400

  # ANALYSED LOG ----
  # load the analysed log
  log_file  <- read.csv(log_path) # %>%print()

  # extract times for dumps
  log_time  <- log_file %>%
    group_by(cond, dump) %>% 
    summarise(time=last(time)) %>% #tail()%>%print()
    rename(seq_time = time, time = dump) # %>%print()
 
  # COVERAGE ----
  # load the coverage data
  cov <- read.csv(input_cov) %>%
    left_join(log_time ) %>%
    arrange(cond, time)%>%print()

  
  if ("bc" %in% colnames(cov)){
    bc_names  <- cov %>%
      filter(bc != "0.0")%>%
      arrange(bc)%>%
      pull(bc)%>%
      unique()
    
    # reorder factor levels
    cov$bc <- factor(cov$bc, levels=bc_names)
  }

  if ("bc" %in% colnames(cov)){
    cov <- cov %>% 
        group_by(cond,bc)
  } else {
    cov <- cov %>% 
        group_by(cond)
  }
    
  # get the ordering of the OTUs
  otu_order <- cov %>%
    filter(cond == 'control', time == max(time)) %>%
    arrange(desc(mean_coverage)) %>%
    pull(otu)%>%
    unique()
  # reorder factor levels
  cov$otu <- factor(cov$otu, levels=otu_order)



  # convert time units to minutes
  cov <- cov %>%
    mutate(seq_time = if_else(is.na(seq_time), 0, seq_time/seq_speed/pores/60),
            evenness = if_else(seq_time==0, NA, evenness))
  


  # UNBLOCK ----
  # load the unblocking data and reorder the factor levels
  unb <- read.csv(input_unb) %>% 
    left_join(log_time) %>% 
    mutate(orig_time = seq_time) %>%
    mutate(seq_time = if_else(is.na(seq_time), 0, seq_time))%>%
    arrange(cond, time)%>%print()
  unb$otu <- factor(unb$otu, levels=otu_order)

  # convert time units to minutes
  unb <- unb %>%
    mutate(seq_time = seq_time/seq_speed/pores/60)

  if ("bc" %in% colnames(cov)){
    unb <- unb %>% 
      group_by(cond, otu)%>%
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


  # SANITY CHECKS ----
  # # Calculate coverage manually from all bases sequenced, not from mapping pileup, as a means of checking results
  # unb %>%
  #   mutate(manual_cov = base_total/genome_size)%>%
  #   select(cond, time, manual_cov)%>%
  #   left_join(., select(cov, cond, time, mean_coverage))

  # # Quantify time
  # alpha = 300
  # rho = 300
  # mu = 400

  # unb %>%
  #   mutate(alphas = total * alpha, rhos = unb * rho) %>% #print(width=Inf)
  #   group_by(cond, time, bc) %>%
  #   summarise(alphas = sum(alphas), rhos = sum(rhos), time_spent_seq = sum(base_total), dump=max(orig_time, na.rm=T)) %>% #print(width=Inf)
  #   mutate(total_time = cumsum(alphas) + cumsum(rhos) + cumsum(time_spent_seq), .before = dump) %>% #print(width=Inf)
  #   mutate(difference = dump - total_time) %>%
  #   print(n=36)

  # unb %>%
  #   mutate(alphas = total * alpha, rhos = unb * rho) %>%
  #   group_by(cond) %>%
  #   summarise(alphas = sum(alphas), rhos = sum(rhos), time_spent_seq = sum(base_total), dump=max(orig_time, na.rm=T)) %>% #print()
  #   mutate(total_time = alphas + rhos + time_spent_seq, .before = dump) %>% #print()
  #   mutate(difference = dump - total_time) %>%
  #   print()

  # Coverage dotplots ----
  # load data
  # filenames_boss <- list.files("temp", pattern="*.boss*.csv", full.names=TRUE)
  # ldf_boss <- lapply(filenames_boss, read.csv)
  # cov_boss <- Reduce(full_join, ldf_boss) %>%
  #   pivot_longer(
  #   !V1,
  #   cols_vary = "fastest",
  #   names_to = "t",
  #   values_to = "cov",
  # )


  # PLOTS ----
  # create plots 
   if ("bc" %in% colnames(cov)){
     col <- "bc"
  } else {
     col <- "otu"
  }

  unb_plot_seq <- ggplot(
      data=unb_cond,
      # mapping=aes(x=seq_time, y=unb_ratio, color=bc)) +
      mapping=aes(x=seq_time, y=unb_ratio)) +
    geom_line(aes(color=.data[[col]])) +
    geom_point(aes(shape=.data[[col]],fill=.data[[col]], size=.data[[col]]), color="white", stroke=0.005, alpha=0.8) +
    scale_color_manual(values=ptol, guide = "none") +
    scale_fill_manual(values=ptol, guide = "none") +
    # scale_linewidth_manual(values=c(1.4,1,0.6))+
    scale_size_manual(values=c(1.8,1.3,1.1), guide = "none")+
    scale_shape_manual(values=c(21,22,23), guide = "none")+
    ylab("cum. rejection rate") +
    xlim(0, max(unb_cond$seq_time))+
    # facet_wrap(~otu, scales="free_y", nrow = nrow) +
    ylim(0,1)+
    theme(legend.position = "none")
  # unb_plot
 ggsave(output, unb_plot_seq, w=4, h=7/3)

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
  # unb_plot
 ggsave(output, unb_plot_frag, w=4, h=7/3)

  unb_noncum_plot_seq <- ggplot(
      data=unb_cond,
      # mapping=aes(x=seq_time, y=unb_ratio, color=bc)) +
      mapping=aes(x=seq_time, y=unb_ratio_noncum, color=.data[[col]])) +
    geom_line() + #linewidth=1,alpha=0.6
    geom_point(aes(shape=.data[[col]],fill=.data[[col]], size=.data[[col]]), color="white", stroke=0.005, alpha=0.8) +
    # ylim(0,1)+
    xlim(0, max(unb_cond$seq_time))+
    scale_color_manual(values=ptol) +
    scale_fill_manual(values=ptol) +
    scale_size_manual(values=c(2,1.4,1))+
    scale_shape_manual(values=c(21,22,23))+
    facet_wrap(~otu, scales="free_y", nrow = nrow) +
    ylab("rejection rate") +
    theme(legend.position = "none")
  # unb_plot
 ggsave(output, unb_noncum_plot_seq, w=4, h=7/3)

  nreads_seq <- ggplot(
      data=unb,
      # mapping=aes(x=seq_time, y=total, linetype=cond, colour=bc)) +
            mapping=aes(x=seq_time, y=total, linetype=cond, linewidth=.data[[col]], colour=.data[[col]])) +

    geom_line() +
    facet_wrap(~otu, scales="free_y", nrow = nrow) +
    xlim(0, max(unb_cond$seq_time))+
    scale_color_manual(values=ptol)+#, guide = "none") +
    scale_linewidth_manual(values=c(1,0.8,0.6))+
    ylab("# fragments")
  # nreads
  ggsave(output, nreads_seq, w=4, h=7/3)
 

  nreads_noncum_seq <- ggplot(
      data=mutate(unb, total_rate = if_else(time == 0, NA, total_rate)),
            mapping=aes(x=seq_time, y=total_rate, linetype=cond, colour=.data[[col]])) +
    geom_line()+
    facet_wrap(~otu, scales="free_y", nrow = nrow) +
    xlim(0, max(unb_cond$seq_time))+
    scale_color_manual(values=ptol) +
    ylab("# fragments/min")+
    theme(legend.position = "none")
  # nreads
  ggsave(output, nreads_noncum_seq, w=4, h=7/3)

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
  ggsave(output, meanc_seq, w=4, h=7/3)


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
  ggsave(output, lowc_seq, w=4, h=7/3)

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
  ggsave(output, evn_seq, w=4, h=7/3)



  time_frag_seq <- ggplot(
      data=unb,
            mapping=aes(x=total, y=seq_time, linetype=cond, colour=.data[[col]], linewidth=.data[[col]])) +
    geom_line() +
    scale_linewidth_manual(values=c(1,0.8,0.6))+
    facet_wrap(~otu, scales="free_y", nrow = nrow) +
    scale_color_manual(values=ptol)+#, guide = "none") +
    ylab("seq. time (minutes)")
  # nreads
  ggsave(output, time_frag_seq, w=4, h=7/3)
  
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
      plot.background = element_rect(fill = "white", color = NA)
    )

  layout_seqt[[2]] <- layout_seqt[[2]] + xlab("# fragments")
  ggsave(output, layout_seqt, w=8, h=7)

}


ptol <- c("#CC6677","#332288","#DDCC77","#117733","#88CCEE","#882255","#44AA99","#999933","#AA4499")
nrow <- 1
visualise_simulation(ptol, nrow)



