#!/usr/bin/env Rscript

library(ggplot2)
library(argparse)
library(dplyr)
library(patchwork)
library(stringr)
library(tidyr)
library(readr)

theme_set(theme_minimal())


get_arguments <- function() {
  parser <- argparse::ArgumentParser()
  parser$add_argument('--input_cov', required = TRUE)
  parser$add_argument('--input_unb', required = TRUE)
  parser$add_argument('--dump_time', required = TRUE)
  parser$add_argument('--pores', required = TRUE)
  parser$add_argument('--analysed_log', required = TRUE)
  parser$add_argument('--output', required = TRUE)
  parser$add_argument('--output_vcf', required = FALSE, default=NULL)
  parser$add_argument('--vcf_summary', required = FALSE, default = NULL, nargs='+')
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

  vcf_summary <- args$vcf_summary
  vcf_output <- args$output_vcf

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
    cov_null_vals <- cov %>% 
      group_by(cond, otu, bc) %>%
      filter(time == 1)%>%
      mutate(time = 0, mean_coverage = 0, low_coverage_prop = 1, evenness = NA, seq_time = 0)

    cov  <- cov %>%
      bind_rows(cov_null_vals)

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

    # get the ordering of the OTUs
    otu_order <- cov %>%
      group_by(cond, bc, otu)%>%
      filter(cond == 'control', time == max(time)) %>% 
      summarise(mean_coverage = max(mean_coverage))%>%
      arrange(bc, desc(mean_coverage)) %>%
      filter(bc == slice_head(ungroup(.), n=1)$bc)%>%
      pull(otu)%>%
      unique()
  } else {
      cov_null_vals <- cov %>% 
        group_by(cond, otu) %>%
        filter(time == 1)%>%
        mutate(time = 0, mean_coverage = 0, low_coverage_prop = 1, evenness = NA, seq_time = 0)

      cov  <- cov %>%
        bind_rows(cov_null_vals)%>% 
        group_by(cond)
      
      # get the ordering of the OTUs
      otu_order <- cov %>%
        filter(cond == 'control', time == max(time)) %>% 
        arrange(desc(mean_coverage)) %>%
        pull(otu)%>%
        unique()
  }
    
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
    unb_null_vals <- unb %>% 
      group_by(cond, otu, bc) %>%
      filter(time == 1)%>%
      mutate(time = 0, total = 0, base_total = 0, unb= 0, unb_ratio = 0, seq_time = 0, orig_time = 0)

    unb  <- unb %>%
      bind_rows(unb_null_vals)%>%
      mutate(orig_time = if_else(is.na(orig_time), 0, orig_time),
      unb_ratio = if_else(is.na(unb_ratio), 0, unb_ratio)) %>%
      arrange(cond,time,bc)

    # reorder factor levels
    unb$bc <- factor(unb$bc, levels=bc_names)

    unb <- unb %>% 
      group_by(cond, bc, otu)
  }
  else{
    unb_null_vals <- unb %>% 
      group_by(cond, otu) %>%
      filter(time == 1)%>%
      mutate(time = 0, total = 0, base_total = 0, unb= 0, unb_ratio = 0, seq_time = 0, orig_time = 0)

    unb  <- unb %>%
      bind_rows(unb_null_vals)%>%
      mutate(orig_time = if_else(is.na(orig_time), 0, orig_time),
      unb_ratio = if_else(is.na(unb_ratio), 0, unb_ratio)) %>%
      arrange(cond,time,bc)%>%
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


  # PLOTS ----
  # create plots 
   if ("bc" %in% colnames(cov)){
     col <- "otu"
     wrap <- "bc"
  } else {
     col <- "otu"
  }

line_thickness <- 0.5
alpha_val <- 0.3

  unb_plot_seq <- ggplot(
      data=unb_cond,
      mapping=aes(x=seq_time, y=unb_ratio), linewidth=line_thickness, alpha=alpha_val) +
    geom_line(aes(color=.data[[col]])) +
    # geom_point(aes(shape=.data[[col]],fill=.data[[col]], size=.data[[col]]), color="white", stroke=0.005, alpha=0.8) +
    scale_color_manual(values=ptol, guide = "none") +
    ylab("cum. rejection rate") +
    xlim(0, max(unb_cond$seq_time))+
    theme(legend.position = "none")+
    ylim(0,1)

  if ("bc" %in% colnames(cov)){
     unb_plot_seq <- unb_plot_seq +
      facet_wrap(~bc, scales="free_y", nrow = nrow)
  }

  # unb_plot
 ggsave(output, unb_plot_seq, w=4, h=7/3)

  unb_plot_frag <- ggplot(
      data=unb_cond,
      mapping=aes(x=total, y=unb_ratio),linewidth=line_thickness, alpha=alpha_val) +
    geom_line(aes(color=.data[[col]])) +
    # geom_point(aes(shape=.data[[col]],fill=.data[[col]], size=.data[[col]]), color="white", stroke=0.005, alpha=0.8) +
    scale_color_manual(values=ptol, guide = "none") +
    ylab("cum. rejection rate") +
    xlim(0, max(unb_cond$total))+
    ylim(0,1)+
    theme(legend.position = "none")

  if ("bc" %in% colnames(cov)){
     unb_plot_frag <- unb_plot_frag +
      facet_wrap(~bc, scales="free_y", nrow = nrow)
  }
  # unb_plot
 ggsave(output, unb_plot_frag, w=4, h=7/3)

  nreads_seq <- ggplot(
      data=unb,
      mapping=aes(x=seq_time, y=total, linetype=cond, colour=.data[[col]]),linewidth=line_thickness, alpha=alpha_val) +
    geom_line() +
    xlim(0, max(unb_cond$seq_time))+
    scale_color_manual(values=ptol)+#, guide = "none") +
    ylab("# fragments")
  # nreads
  if ("bc" %in% colnames(cov)){
     nreads_seq <- nreads_seq +
      facet_wrap(~bc, scales="free_y", nrow = nrow)
  }
  ggsave(output, nreads_seq, w=4, h=7/3)


  meanc_seq <- cov %>% 
    filter(seq_time <= max(unb_cond$seq_time))%>%
    ggplot( mapping=aes(x=seq_time, y=mean_coverage, linetype=cond, colour=.data[[col]]),linewidth=line_thickness, alpha=alpha_val) +
    geom_line() +
    # facet_wrap(~otu, scales="free_y", nrow = nrow) +
    xlim(0, max(unb_cond$seq_time))+
    # ylim(0, max(pull(slice_min(filter(cov, seq_time > max(unb_cond$seq_time)), order_by=time), mean_coverage)))+
    scale_color_manual(values=ptol)+ #, guide = "none") +
    # scale_linewidth_manual(values=c(1,0.8,0.6))+
    # scale_alpha_manual(values=c(1,0.9,0.8), guide = "none")+
    ylab("mean coverage")+
    theme(legend.position = "none")
  if ("bc" %in% colnames(cov)){
     meanc_seq <- meanc_seq +
      facet_wrap(~bc, scales="free_y", nrow = nrow)
  }
  # meanc
  ggsave(output, meanc_seq, w=4, h=7/3)


  lowc_seq <- ggplot(
      data=cov,
            mapping=aes(x=seq_time, y=low_coverage_prop, linetype=cond, colour=.data[[col]]),linewidth=line_thickness, alpha=alpha_val) +
    geom_line() +
    # facet_wrap(~otu, scales="free_y", nrow = nrow) +
    scale_color_manual(values=ptol) +
    xlim(0, max(unb_cond$seq_time))+
    # scale_linewidth_manual(values=c(1,0.8,0.6))+
    ylab("prop. sites at <5x")+
    theme(legend.position = "none")
  # lowc
  if ("bc" %in% colnames(cov)){
     lowc_seq <- lowc_seq +
      facet_wrap(~bc, scales="free_y", nrow = nrow)
  }
  ggsave(output, lowc_seq, w=4, h=7/3)

  evn_seq <- ggplot(
      data=cov,
            mapping=aes(x=seq_time, y=evenness, linetype=cond, colour=.data[[col]]),linewidth=line_thickness, alpha=alpha_val) +
    geom_line() +
    # facet_wrap(~otu, scales="free_y", nrow = nrow) +
    scale_color_manual(values=ptol) +
    xlim(0, max(unb_cond$seq_time))+
    # scale_linewidth_manual(values=c(1,0.8,0.6))+
    # scale_alpha_manual(values=c(0.6,1), guide = "none")+
    ylab("evenness")+
    theme(legend.position = "none")

  if ("bc" %in% colnames(cov)){
     evn_seq <- evn_seq +
      facet_wrap(~bc, scales="free_y", nrow = nrow)
  }
  # evenness of coverage
  ggsave(output, evn_seq, w=4, h=7/3)

  
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
  
  
  # Variant call summary plots if relevant
  if (vcf_summary != NULL){
    # load vcf summary file
    processed_hap <- read_csv(vcf_summary) %>%
      rename("cond" = "exp", "time" = "time_point", "bc"="barcode")%>%
      left_join(cov)%>%
      mutate(seq_time = seq_time/24)

    precision_seqt <- processed_hap %>%
      filter(Filter=="PASS", Type=="SNP")%>%
      ggplot()+
      geom_line(aes(x=seq_time, y=METRIC.Precision, group=cond, linetype=cond, colour=otu))+
      ylim(0,1)+
      facet_grid(rows=vars(otu))

    recall_seqt <- processed_hap %>%
      filter(Filter=="PASS", Type=="SNP")%>%
      ggplot()+
      geom_line(aes(x=seq_time, y=METRIC.Recall, group=cond, linetype=cond, colour=otu))+
      ylim(0,1)+
      facet_grid(rows=vars(otu))

    precision_mcov <- processed_hap %>%
      filter(Filter=="PASS", Type=="SNP")%>%
      ggplot()+
      geom_line(aes(x=mean_coverage, y=METRIC.Precision, group=cond, linetype=cond, colour=otu))+
      ylim(0,1)+
      facet_grid(rows=vars(otu))

    recall_mcov <- processed_hap %>%
      filter(Filter=="PASS", Type=="SNP")%>%
      ggplot()+
      geom_line(aes(x=mean_coverage, y=METRIC.Recall, group=cond, linetype=cond, colour=otu))+
      ylim(0,1)+
      facet_grid(rows=vars(otu))


    layout_vcf <- wrap_plots(list(recall_seqt, recall_mcov,precision_seqt,precision_mcov), axes='keep', axis_titles="keep")+
      plot_annotation(tag_levels = "A") +
      plot_layout(
        nrow = 1,
        ncol = 4,
        guides="collect"
        ) &
      xlab("seq. time (pseudohours)") &
      guides(color = NULL)&
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.background = element_rect(fill = "white", color = NA)
      )

    layout_vcf[[2]] <- layout_vcf[[2]] + xlab("mean coverage")
    layout_vcf[[4]] <- layout_vcf[[4]] + xlab("mean coverage")

    ggsave(vcf_output, layout_vcf, w=9, h=6)
  }


}


ptol <- c("#CC6677","#332288","#DDCC77","#117733","#88CCEE","#882255","#44AA99","#999933","#AA4499")
nrow <- 1
visualise_simulation(ptol, nrow)


