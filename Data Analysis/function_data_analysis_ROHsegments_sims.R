
require(ggplot2)
require(data.table)
require(ggpubr)
require(vroom)
require(dplyr)
require(OneR)

script_path = '/Users/gabrielemariasgarlata/Documents/Projects/Ongoing_projects/Response2Orkin_et_al/bin/IBDsegments_Panmictic'
source(paste0(script_path,'/Functions_Panmictic_PopSizeChange.R'))
source(paste0(script_path,'/Functions_Panmictic_ConstantNe.R'))


GetSimulatedROHsegments_BinnedByAge_Across_Scenarios<-function(pre_pathdata,n_sims,n_bin_age,binwidth_val,selectedPAIR,G,rbp,SampleSize,ConstantNe,TimePopChange,AncientNe_expansion,RecentNe_expansion,AncientNe_bottleneck,RecentNe_bottleneck){
  
  Lcenti = 0.01/rbp
  G_M = G*rbp
  
  constant_scenario = paste0('Panmictic_Ne',ConstantNe,'_SampleSize',SampleSize,'_rbp',rbp,'_GenomeLength',as.character(G))
  pre_inputname_constant = paste0('IBDsegments_',constant_scenario)
  path_SimulatedIBDsegm_constant = paste0(pre_pathdata,'/',constant_scenario)
  
  intervalDF_constant = MakeTheoreticalSimulatedDF_OneIndAcrossSimsMinMaxInterval(pre_inputname = pre_inputname_constant, selectedPAIR = selectedPAIR, Lcenti = Lcenti,
                                                                                  binwidth=binwidth_val, N = ConstantNe, n_sims = n_sims, path_SimulatedIBDsegm = path_SimulatedIBDsegm_constant)
  
  Gtot_constant = intervalDF_constant$`Genome Length`
  resDF_constant = intervalDF_constant$ExpSimDF
  rawDataDF_constant = intervalDF_constant$`Raw Data`
  
  resDF_constant$x_cM = resDF_constant$x * 100
  
  ageDF_SimulatedIBDsegm_Constant = GetAgeSimulated_IBDsegments(rawDataDF = rawDataDF_constant, resDF = resDF_constant)
  ageDF_SimulatedIBDsegm_Constant$type = 'Constant'
  ##END: Constant Pop Size##
  
  
  ##START: Bottleneck Pop Size##
  bottleneck_scenario = paste0('Panmictic_ancientNe',AncientNe_bottleneck,'_PopSizeChangeTime',TimePopChange,'_recentNe',RecentNe_bottleneck,'_SampleSize',SampleSize,'_rbp',rbp,'_GenomeLength',as.character(G))
  pre_inputname_bottleneck = paste0('IBDsegments_',bottleneck_scenario)
  path_SimulatedIBDsegm_bottleneck = paste0(pre_pathdata,'/',bottleneck_scenario)
  
  
  
  intervalDF_bottleneck = MakeTheoreticalSimulatedDF_OneIndAcrossSimsMinMaxInterval_PopChange(pre_inputname = pre_inputname_bottleneck, selectedPAIR = selectedPAIR, Lcenti = Lcenti,
                                                                                              binwidth=binwidth_val, t_change = TimePopChange, recentNe = RecentNe_bottleneck, 
                                                                                              ancientNe = AncientNe_bottleneck, n_sims = n_sims, path_SimulatedIBDsegm = path_SimulatedIBDsegm_bottleneck)
  
  
  Gtot_bottleneck = intervalDF_bottleneck$`Genome Length`
  resDF_bottleneck = intervalDF_bottleneck$ExpSimDF
  rawDataDF_bottleneck = intervalDF_bottleneck$`Raw Data`
  resDF_bottleneck$x_cM = resDF_bottleneck$x * 100
  
  ageDF_SimulatedIBDsegm_Bottleneck = GetAgeSimulated_IBDsegments(rawDataDF = rawDataDF_bottleneck, resDF = resDF_bottleneck)
  ageDF_SimulatedIBDsegm_Bottleneck$type = 'Bottleneck'
  ##END: Bottleneck Pop Size##
  
  
  ##START: Expansion Pop Size##
  
  expansion_scenario = paste0('Panmictic_ancientNe',AncientNe_expansion,'_PopSizeChangeTime',TimePopChange,'_recentNe',RecentNe_expansion,'_SampleSize',SampleSize,'_rbp',rbp,'_GenomeLength',as.character(G))
  pre_inputname_expansion = paste0('IBDsegments_',expansion_scenario)
  path_SimulatedIBDsegm_expansion = paste0(pre_pathdata,'/',expansion_scenario)
  
  
  
  intervalDF_expansion = MakeTheoreticalSimulatedDF_OneIndAcrossSimsMinMaxInterval_PopChange(pre_inputname = pre_inputname_expansion, selectedPAIR = selectedPAIR, Lcenti = Lcenti,
                                                                                             binwidth=binwidth_val, t_change = TimePopChange, recentNe = RecentNe_expansion, 
                                                                                             ancientNe = AncientNe_expansion, n_sims = n_sims, path_SimulatedIBDsegm = path_SimulatedIBDsegm_expansion)
  
  
  Gtot_expansion = intervalDF_expansion$`Genome Length`
  resDF_expansion = intervalDF_expansion$ExpSimDF
  rawDataDF_expansion = intervalDF_expansion$`Raw Data`
  resDF_expansion$x_cM = resDF_expansion$x * 100
  
  ageDF_SimulatedIBDsegm_Expansion = GetAgeSimulated_IBDsegments(rawDataDF = rawDataDF_expansion, resDF = resDF_expansion)
  ageDF_SimulatedIBDsegm_Expansion$type = 'Expansion'
  ##END: Expansion Pop Size##
  
  Tot_ageDF_SimulatedIBDsegm = NULL
  Tot_ageDF_SimulatedIBDsegm = rbind(Tot_ageDF_SimulatedIBDsegm,ageDF_SimulatedIBDsegm_Constant)
  Tot_ageDF_SimulatedIBDsegm = rbind(Tot_ageDF_SimulatedIBDsegm,ageDF_SimulatedIBDsegm_Bottleneck)
  Tot_ageDF_SimulatedIBDsegm = rbind(Tot_ageDF_SimulatedIBDsegm,ageDF_SimulatedIBDsegm_Expansion)
  
  
  age_binnedDF = Tot_ageDF_SimulatedIBDsegm %>%
    reframe(binID = bin(IBDage, nbins = n_bin_age, method="length"), IBDage = IBDage, x_cM = x_cM, xmin = xmin, xmax = xmax, type = type)
  
  finalPlotAgeDF = age_binnedDF %>%
    group_by(x_cM, binID,type) %>%
    summarise(meanTage = mean(IBDage), count = length(IBDage))
  
  tmpBinLabel = as.data.frame(do.call(rbind,strsplit(as.character(finalPlotAgeDF$binID),',')))
  
  minBinAge = round(as.numeric(gsub('\\(','',tmpBinLabel$V1)))
  minBinAge[which(minBinAge<0)] = 0
  maxBinAge = round(as.numeric(gsub('\\]','',tmpBinLabel$V2)))
  
  finalPlotAgeDF$time_min = minBinAge
  finalPlotAgeDF$time_max = maxBinAge
  
  finalPlotAgeDF = finalPlotAgeDF[order(finalPlotAgeDF$time_min),]
  finalPlotAgeDF$binLabel = paste0(finalPlotAgeDF$time_min,'-',finalPlotAgeDF$time_max)
  
  finalPlotAgeDF$binLabel = factor(finalPlotAgeDF$binLabel, levels = unique(finalPlotAgeDF$binLabel))
  
  finalPlotAgeDF_constant = subset(finalPlotAgeDF,type=='Constant')
  finalPlotAgeDF_expansion = subset(finalPlotAgeDF,type=='Expansion')
  finalPlotAgeDF_bottleneck = subset(finalPlotAgeDF,type=='Bottleneck')
  
  maxval_ylimDF = finalPlotAgeDF %>% 
    group_by(type, x_cM) %>%
    summarise(tot_count = sum(count))
  
  maxval_ylim = max(maxval_ylimDF$tot_count)
  
  
  plotSimDataByAgeAndBlockBin_Constant = ggplot(finalPlotAgeDF_constant,aes(x=x_cM,y=count,group = binLabel,colour=binLabel,fill=binLabel))+
    geom_bar(position = "stack",stat='identity')+
    xlab('cM')+
    ylab('#IBD segments')+
    scale_fill_discrete(name = 'Bin Age')+
    scale_color_discrete(name = 'Bin Age')+
    xlim(0,0.25)+
    ylim(0,maxval_ylim)+
    ggtitle('Constant')+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))
  
  plotSimDataByAgeAndBlockBin_Bottleneck = ggplot(finalPlotAgeDF_bottleneck,aes(x=x_cM,y=count,group = binLabel,colour=binLabel,fill=binLabel))+
    geom_bar(position = "stack",stat='identity')+
    xlab('cM')+
    ylab('#IBD segments')+
    scale_fill_discrete(name = 'Bin Age')+
    scale_color_discrete(name = 'Bin Age')+
    xlim(0,0.25)+
    ylim(0,maxval_ylim)+
    ggtitle('Bottleneck')+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))
  
  plotSimDataByAgeAndBlockBin_Expansion = ggplot(finalPlotAgeDF_expansion,aes(x=x_cM,y=count,group = binLabel,colour=binLabel,fill=binLabel))+
    geom_bar(position = "stack",stat='identity')+
    xlab('cM')+
    ylab('#IBD segments')+
    scale_fill_discrete(name = 'Bin Age')+
    scale_color_discrete(name = 'Bin Age')+
    xlim(0,0.25)+
    ylim(0,maxval_ylim)+
    ggtitle('Expansion')+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  finalPlot = ggarrange(print(plotSimDataByAgeAndBlockBin_Constant), print(plotSimDataByAgeAndBlockBin_Bottleneck), 
                        print(plotSimDataByAgeAndBlockBin_Expansion), ncol = 3, nrow = 1, common.legend = TRUE, legend = 'right')
  
  return(finalPlot)
}


GetChronROH_trueROHsegments_vs_trueFroh<-function(pre_path_data,G_bp,n_sims,TimePopChange,SamplingTime_list,scenario_list, labels_list){
  
  
  index = 1
  
  for(scenario_name in scenario_list){
    
    label_plot = labels_list[index]
    
    path_chrono_roh_data = paste0(pre_path_data,'/',scenario_name)
    
    final_bcftools_roh = NULL
    final_chrono_roh = NULL
    
    for(simID in 1:n_sims){
      
      df_chrono = fread(file = paste0(path_chrono_roh_data,'/',scenario_name,'_time0_chronoRoh_sim',simID,'.csv'))
      
      df_chrono$type = factor(df_chrono$type, levels = unique(df_chrono$type))
      
      path_roh_data = paste0(pre_path_data,'/',scenario_name,'/sim',simID)
      
      final_roh_sim = NULL
      
      for(SamplingTime in SamplingTime_list){
        
        input_filename = paste0(scenario_name,'_Time',SamplingTime,'_roh.txt')
        
        input_file = paste0(path_roh_data,'/',input_filename)
        
        df_ROH = read.table(input_file)
        
        names(df_ROH) = c('RG','sample','chromosome',	'start', 'end', 'length_bp', 'n_markers',	'quality')
        
        totROH_per_ind = df_ROH %>%
          group_by(sample) %>%
          summarise(l_roh = sum(length_bp))
        
        F_roh_val = mean(totROH_per_ind$l_roh)/G_bp
        
        tmp_roh = data.frame(time = SamplingTime, F_roh = F_roh_val)
        
        final_roh_sim = rbind(final_roh_sim,tmp_roh)
      }
      
      final_roh_sim$simID = simID
      final_bcftools_roh = rbind(final_bcftools_roh,final_roh_sim)
      
      df_chrono$simID = simID
      df_chrono$combo = paste0(df_chrono$type,'_',df_chrono$ind,'_',df_chrono$simID)
      
      final_chrono_roh = rbind(final_chrono_roh,df_chrono)
      
    }
    
    
    final_bcftools_roh = data.table(age = final_bcftools_roh$time, F_roh = final_bcftools_roh$F_roh, type = 'sampling_time',
                                    ind = 'avg', simID = final_bcftools_roh$simID)
    
    final_bcftools_roh$combo = paste0(final_bcftools_roh$ind,'_',final_bcftools_roh$simID)
    
    final_chrono_roh = rbind(final_chrono_roh,final_bcftools_roh)
    
    list_combo = unique(final_chrono_roh$combo)
    
    list_combo_no_avg = list_combo[-grep('avg_',list_combo)]
    list_combo_avg = list_combo[grep('avg_',list_combo)]
    
    
    final_chrono_roh$combo = factor(final_chrono_roh$combo, levels = c(list_combo_no_avg,list_combo_avg))
    
    final_chrono_roh_estIBDage = final_chrono_roh[!final_chrono_roh$type=='IBDage',]
    
    plotChronoRoh_estIBDage = ggplot()+
      geom_line(data = final_chrono_roh_estIBDage, aes(x = age, y = F_roh, group = combo, linetype = type, colour = type))+
      #geom_point(data = final_chrono_roh_estIBDage, aes(x = age, y = F_roh, group = combo, linetype = type, colour = type,alpha = type), size =1)+
      xlab('Generations')+ylab(expression(F["ROH"]))+
      scale_linetype_manual(name = '',values = c(1, 1), labels = c('estIBDage' = '   Estimated Age\n(as in Orkin et al., 2025)', 'sampling_time' = expression(True~F["ROH"])))+
      scale_colour_manual(name = '', values = c("#E69F00","#0072B2"),labels = c('estIBDage' = '   Estimated Age\n(as in Orkin et al., 2025)', 'sampling_time' = expression(True~F["ROH"])))+
      #scale_size_manual(name = '', values = c(.2,.5), labels = c('estIBDage' = '   Estimated Age\n(as in Orkin et al., 2025)', 'sampling_time' = expression(True~F["ROH"])),guide = 'none')+
      #scale_alpha_manual(name = '', values = c(.3,1), labels = c('estIBDage' = '   Estimated Age\n(as in Orkin et al., 2025)', 'sampling_time' = expression(True~F["ROH"])),guide = 'none')+
      xlim(0,1000)+
      theme_bw()+
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5)) 
    

    if(label_plot=='Constant'){
      
      plotChronoRoh_estIBDage = plotChronoRoh_estIBDage+
        ggtitle('Constant')+theme(axis.line = element_line(colour = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  panel.border = element_blank(),
                                  panel.background = element_blank(),
                                  plot.title = element_text(hjust = 0.5))
    }
    
    if(label_plot=='Expansion'){
      plotChronoRoh_estIBDage = plotChronoRoh_estIBDage + geom_vline(xintercept = TimePopChange, linetype = 'dashed')+
        ggtitle('Expansion')+theme(axis.line = element_line(colour = "black"),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   panel.border = element_blank(),
                                   panel.background = element_blank(),
                                   plot.title = element_text(hjust = 0.5))
    }
    
    
    if(label_plot=='Bottleneck'){

      plotChronoRoh_estIBDage = plotChronoRoh_estIBDage + geom_vline(xintercept = TimePopChange, linetype = 'dashed')+
        ggtitle('Bottleneck')+theme(axis.line = element_line(colour = "black"),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    panel.border = element_blank(),
                                    panel.background = element_blank(),
                                    plot.title = element_text(hjust = 0.5))
    }
    
    assign(paste0('chronoEstAge_plot', index),plotChronoRoh_estIBDage)
    
    index = index + 1
    
  }
  
  finalPlot_EstAge = ggarrange(print(chronoEstAge_plot1), print(chronoEstAge_plot2), print(chronoEstAge_plot3), ncol = 3, nrow = 1, common.legend = TRUE, legend = 'right')
  

  return(finalPlot_EstAge)
}


GetChronROH_inferredROHsegments_vs_trueFroh<-function(pre_path_data,G_bp,n_sims,TimePopChange,SamplingTime_list,scenario_list, labels_list, bin_size, standardize_Froh){
  
  
  index = 1
  
  for(scenario_name in scenario_list){
    
    label_plot = labels_list[index]
    
    path_chrono_roh_data = paste0(pre_path_data,'/',scenario_name)
    
    final_bcftools_roh = NULL
    final_chrono_roh = NULL
    
    for(simID in 1:n_sims){
      
      df_chrono = fread(file = paste0(path_chrono_roh_data,'/',scenario_name,'_time0_chronoRoh_bcftools_sim',simID,'_bin_interval',bin_size,'.csv'))
      
      df_chrono$type = factor(df_chrono$type, levels = unique(df_chrono$type))
      
      path_roh_data = paste0(pre_path_data,'/',scenario_name,'/sim',simID)
      
      final_roh_sim = NULL
      
      for(SamplingTime in SamplingTime_list){
        
        input_filename = paste0(scenario_name,'_Time',SamplingTime,'_roh.txt')
        
        input_file = paste0(path_roh_data,'/',input_filename)
        
        df_ROH = read.table(input_file)
        
        names(df_ROH) = c('RG','sample','chromosome',	'start', 'end', 'length_bp', 'n_markers',	'quality')
        
        totROH_per_ind = df_ROH %>%
          group_by(sample) %>%
          summarise(l_roh = sum(length_bp))
        
        F_roh_val = mean(totROH_per_ind$l_roh)/G_bp
        
        tmp_roh = data.frame(time = SamplingTime, F_roh = F_roh_val)
        
        final_roh_sim = rbind(final_roh_sim,tmp_roh)
      }
      
      final_roh_sim$simID = simID
      final_bcftools_roh = rbind(final_bcftools_roh,final_roh_sim)
      
      df_chrono$simID = simID
      df_chrono$combo = paste0(df_chrono$type,'_',df_chrono$ind,'_',df_chrono$simID)
      
      final_chrono_roh = rbind(final_chrono_roh,df_chrono)
      
    }
    
    
    final_bcftools_roh = data.table(age = final_bcftools_roh$time, F_roh = final_bcftools_roh$F_roh, type = 'sampling_time',
                                    ind = 'avg', simID = final_bcftools_roh$simID)
    
    final_bcftools_roh$combo = paste0(final_bcftools_roh$ind,'_',final_bcftools_roh$simID)
    
    final_chrono_roh = rbind(final_chrono_roh,final_bcftools_roh)
    
    list_combo = unique(final_chrono_roh$combo)
    
    list_combo_no_avg = list_combo[-grep('avg_',list_combo)]
    list_combo_avg = list_combo[grep('avg_',list_combo)]
    
    
    final_chrono_roh$combo = factor(final_chrono_roh$combo, levels = c(list_combo_no_avg,list_combo_avg))
    
    if(standardize_Froh==TRUE){
      final_chrono_roh_estROHage = subset(final_chrono_roh,type=='estIBDage_bcftools')
      final_chrono_roh_estROHage$F_roh =  scale(final_chrono_roh_estROHage$F_roh)
      
      final_chrono_roh_trueFroh = subset(final_chrono_roh,type=='sampling_time')
      final_chrono_roh_trueFroh$F_roh =  scale(final_chrono_roh_trueFroh$F_roh)
      
      final_chrono_roh = rbind(final_chrono_roh_estROHage,final_chrono_roh_trueFroh)
    }else{}
    
    
    plotChronoRoh_estIBDage = ggplot()+
      geom_line(data = final_chrono_roh, aes(x = age, y = F_roh, group = combo, linetype = type, colour = type))+
      #geom_point(data = final_chrono_roh_estIBDage, aes(x = age, y = F_roh, group = combo, linetype = type, colour = type,alpha = type), size =1)+
      xlab('Generations')+
      scale_linetype_manual(name = '',values = c(1, 1), labels = c('estIBDage_bcftools' = '   Estimated Age\n(as in Orkin et al., 2025)', 'sampling_time' = expression(True~F["ROH"])))+
      scale_colour_manual(name = '', values = c("#E69F00","#0072B2"),labels = c('estIBDage_bcftools' = '   Estimated Age\n(as in Orkin et al., 2025)', 'sampling_time' = expression(True~F["ROH"])))+
      scale_size_manual(name = '', values = c(.2,.5), labels = c('estIBDage_bcftools' = '   Estimated Age\n(as in Orkin et al., 2025)', 'sampling_time' = expression(True~F["ROH"])),guide = 'none')+
      #scale_alpha_manual(name = '', values = c(.3,1), labels = c('estIBDage_bcftools' = '   Estimated Age\n(as in Orkin et al., 2025)', 'sampling_time' = expression(True~F["ROH"])),guide = 'none')+
      xlim(0,1000)+
      theme_bw()+
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5)) 
    
    
    if(label_plot=='Constant'){
      
      plotChronoRoh_estIBDage = plotChronoRoh_estIBDage+
        ggtitle('Constant')+theme(axis.line = element_line(colour = "black"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  panel.border = element_blank(),
                                  panel.background = element_blank(),
                                  plot.title = element_text(hjust = 0.5))
    }
    
    if(label_plot=='Expansion'){
      
      plotChronoRoh_estIBDage = plotChronoRoh_estIBDage + geom_vline(xintercept = TimePopChange, linetype = 'dashed')+
        ggtitle('Expansion')+theme(axis.line = element_line(colour = "black"),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   panel.border = element_blank(),
                                   panel.background = element_blank(),
                                   plot.title = element_text(hjust = 0.5))
    }
    
    
    if(label_plot=='Bottleneck'){
      
      plotChronoRoh_estIBDage = plotChronoRoh_estIBDage + geom_vline(xintercept = TimePopChange, linetype = 'dashed')+
        ggtitle('Bottleneck')+theme(axis.line = element_line(colour = "black"),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    panel.border = element_blank(),
                                    panel.background = element_blank(),
                                    plot.title = element_text(hjust = 0.5))
    }
    
    if(standardize_Froh==TRUE){
      plotChronoRoh_estIBDage = plotChronoRoh_estIBDage + ylab(expression(standardized~F["ROH"]))
    }else{
      plotChronoRoh_estIBDage = plotChronoRoh_estIBDage + ylab(expression(F["ROH"]))
    }
    
    assign(paste0('chronoEstAge_plot', index),plotChronoRoh_estIBDage)
    
    index = index + 1
    
  }
  
  finalPlot_EstAge = ggarrange(print(chronoEstAge_plot1), print(chronoEstAge_plot2), print(chronoEstAge_plot3), ncol = 3, nrow = 1, common.legend = TRUE, legend = 'right')
  
  return(finalPlot_EstAge)
}

