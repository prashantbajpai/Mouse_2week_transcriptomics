awesomefunction = function(res, title = '', my_xmin = -1.5, my_xmax = 1.5, 
                           lab = res$GeneID, ...){
  up_genes = res[which(res$avg_log2FC >= 0.378 & res$p_val_adj <= 0.05), ]
  down_genes = res[which(res$avg_log2FC <= -0.378 & res$p_val_adj <= 0.05), ]
  #Volcano Plot
  library(RColorBrewer)
  library(pheatmap)
  library(EnhancedVolcano)
  #library(xlsx)
  library(ggpubr)
  library(scales)
  res = res[!is.na(res$p_val_adj),]
  res$p_val_adj[res$p_val_adj == 0] = min(res$p_val_adj[res$p_val_adj != 0])*0.1
  #res$p_val_adj[-log10(res$p_val_adj) > 20] = rescale(res$p_val_adj[-log10(res$p_val_adj) > 20],
  #                                          c(1e-16, 1e-20))
  #res$avg_log2FC[res$avg_log2FC > 10] = NA
  #res$avg_log2FC[res$avg_log2FC < -10] = NA
  keyval = ifelse(res$avg_log2FC >= 0.378 & res$p_val_adj <= 0.05, '#b2182b', 
                  ifelse(res$avg_log2FC <= -0.378 & res$p_val_adj <= 0.05, '#2166ac',
                         'grey70'))
  names(keyval)[keyval == '#b2182b'] <- 'high'
  names(keyval)[keyval == '#2166ac'] <- 'mid'
  names(keyval)[keyval == 'grey70'] <- 'low'
  #my_xmin = min(res$avg_log2FC, na.rm = T)
  #my_xmax = max(res$avg_log2FC, na.rm = T)
  my_ymax = max(-log10(res$p_val_adj), na.rm = T)
  p1 = EnhancedVolcano(res, lab = lab, 
                       x = 'avg_log2FC', y = 'p_val_adj', 
                       xlab = bquote(~Log[2]~ 'fold change'),
                       ylab = bquote(~-Log[10]~ 'adjusted P'),
                       pCutoff = 0.05, FCcutoff = 0.378, pointSize = 1.5,
                       #xlim = c(-10, 10),
                       ylim = c(0, my_ymax*1.1),
                       labSize = 4.0, colAlpha = 0.6, legendPosition = 'none',
                       subtitle = '',
                       title = title,
                       caption = '',
                       legendLabSize = 12, legendIconSize = 4.0,
                       titleLabSize = 24,
                       drawConnectors = T,
                       colCustom = keyval,
                       gridlines.major = F,
                       gridlines.minor = F,
                       max.overlaps = 10)
  up_gene_count = paste(dim(up_genes)[1])
  down_gene_count = paste(dim(down_genes)[1])
  p1 = p1 + annotate('text', label = up_gene_count, x = my_xmax, y = my_ymax*1.1, size = 7, 
                     color = '#b2182b', fontface = 'bold') + 
    annotate('text', label = down_gene_count, x = my_xmin, y = my_ymax*1.1, size = 7, 
             color = '#2166ac', fontface = 'bold') +
    theme(plot.margin = margin(t = 0, r = 10, b = -20, l = 10),
          plot.title = element_text(hjust = 0.5, vjust = 0))
  return(p1)
}

PrepareCellchatdata <- function(seuratObj, group, groupby, species, 
                                groupcol, pvalthresh = 0.05){
  data.input = GetAssayData(seuratObj, assay = "RNA", slot = "data") 
  meta = seuratObj@meta.data
  
  colnames(meta)[colnames(meta) == groupcol] = "Group"
  cell.use = rownames(meta)[meta$Group == group]
  data.input.fil = data.input[,cell.use]
  meta.fil = meta[cell.use,]
  message(paste("filtering cells:", unique(meta.fil$Group)))
  
  cellchat = createCellChat(object = data.input.fil, meta = meta.fil, 
                            group.by = groupby)
  
  message("Setting the ligand-receptor interaction")
  if(species == "human"){
    CellChatDB = CellChatDB.human
  }else if(species == "mouse"){
    CellChatDB = CellChatDB.mouse
    #showDatabaseCategory(CellChatDB)
    CellChatDB[["interaction"]] = CellChatDB[["interaction"]][-1887,]
    CellChatDB[["interaction"]] = CellChatDB[["interaction"]][-1900,]
  }else{
    stop("Incorrect species name, valid entries: mouse, human")
  }
  
  cellchat@DB = CellChatDB
  cellchat = subsetData(cellchat, features = NULL)
  future::plan("multisession", workers = 8)
  cellchat = identifyOverExpressedGenes(cellchat)
  cellchat = identifyOverExpressedInteractions(cellchat)
  if(species == "human"){
    cellchat = projectData(cellchat, PPI.human)
  }else{
    cellchat = projectData(cellchat, PPI.mouse)
  }
  message("Computing communication probability and infering cellular communication network")
  cellchat = computeCommunProb(cellchat, type = 'triMean', population.size = T)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat, thresh = pvalthresh)
  cellchat <- aggregateNet(cellchat, thresh = pvalthresh)
  return(cellchat)
}

pathwayprobsum <- function (object, signaling, signaling.name = NULL, 
                            vertex.weight = 1, vertex.size.max = NULL, thresh = 0.05, 
                            ...) {
  
  if (is.null(vertex.weight)) {
    vertex.weight <- as.numeric(table(object@idents))
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    }
    else {
      vertex.size.max <- 15
    }
  }
  pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, 
                       key = "pathway_name", matching.exact = T, pair.only = T)
  if (is.null(signaling.name)) {
    signaling.name <- signaling
  }
  net <- object@net
  pairLR.use.name <- dimnames(net$prob)[[3]]
  pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
  pairLR <- pairLR[pairLR.name, ]
  prob <- net$prob
  pval <- net$pval
  prob[pval > thresh] <- 0
  if (length(pairLR.name) > 1) {
    pairLR.name.use <- pairLR.name[apply(prob[, , pairLR.name], 
                                         3, sum) != 0]
  }
  else {
    pairLR.name.use <- pairLR.name[sum(prob[, , pairLR.name]) != 
                                     0]
  }
  if (length(pairLR.name.use) == 0) {
    stop(paste0("There is no significant communication of ", 
                signaling.name))
  }
  else {
    pairLR <- pairLR[pairLR.name.use, ]
  }
  nRow <- length(pairLR.name.use)
  prob <- prob[, , pairLR.name.use]
  pval <- pval[, , pairLR.name.use]
  if (length(dim(prob)) == 2) {
    prob <- replicate(1, prob, simplify = "array")
    pval <- replicate(1, pval, simplify = "array")
  }
  
  prob.sum <- apply(prob, c(1, 2), sum)
  
  return(prob.sum)
}

# Function to pseudo customize FeaturePlots
customize_Seurat_FeaturePlot <- function(p, alpha.use = 1, gradient.use = c("yellow", "red"), expression.threshold = NULL, is.log1p.transformed = F) {
  
  #### Main function ####
  main_function <- function(p = p, alpha.use = alpha.use, gradient.use = gradient.use, expression.threshold = expression.threshold, is.log1p.transformed = is.log1p.transformed) {
    
    # Remove cells having an expression below a certain threshold
    if (!is.null(expression.threshold)) {
      if (isTRUE(is.log1p.transformed)) {
        p$data <- p$data[p$data$gene >= expression.threshold,]
      } else {
        p$data <- p$data[p$data$gene >= log1p(expression.threshold),]
      }
    }
    
    # Fill points using the gene expression levels
    p$layers[[1]]$mapping$fill <- p$layers[[1]]$mapping$colour
    
    # Define transparency of points
    p$layers[[1]]$mapping$alpha <- alpha.use
    
    # Change fill and colour gradient values
    p <- p + scale_colour_gradientn(colours = gradient.use, guide = F) +
      scale_fill_gradientn(colours = gradient.use, name = expression(atop(Expression, (log)))) +
      scale_alpha_continuous(range = alpha.use, guide = F)
  }
  
  #### Execution of main function ####
  # Apply main function on all features
  p <- lapply(X = p, alpha.use = alpha.use, gradient.use = gradient.use, 
              expression.threshold = expression.threshold, is.log1p.transformed = is.log1p.transformed,
              FUN = main_function)
  
  # Arrange all plots using cowplot
  # Adapted from Seurat
  # https://github.com/satijalab/seurat/blob/master/R/plotting.R#L1100
  # ncol argument from Josh O'Brien
  # https://stackoverflow.com/questions/10706753/how-do-i-arrange-a-variable-list-of-plots-using-grid-arrange
  cowplot::plot_grid(plotlist = p, ncol = ceiling(sqrt(length(p))))
}

awesomefunction2 = function(res, title = '', my_xmin = -1.5, my_xmax = 1.5, 
                           lab = res$GeneID, xlim = c(-2,2), ...){
  up_genes = res[which(res$avg_log2FC >= 0.378 & res$p_val_adj <= 0.05), ]
  down_genes = res[which(res$avg_log2FC <= -0.378 & res$p_val_adj <= 0.05), ]
  #Volcano Plot
  library(RColorBrewer)
  library(pheatmap)
  library(EnhancedVolcano)
  #library(xlsx)
  library(ggpubr)
  library(scales)
  res = res[!is.na(res$p_val_adj),]
  res$p_val_adj[res$p_val_adj == 0] = min(res$p_val_adj[res$p_val_adj != 0])*0.1
  #res$p_val_adj[-log10(res$p_val_adj) > 20] = rescale(res$p_val_adj[-log10(res$p_val_adj) > 20],
  #                                          c(1e-16, 1e-20))
  #res$avg_log2FC[res$avg_log2FC > 10] = NA
  #res$avg_log2FC[res$avg_log2FC < -10] = NA
  keyval = ifelse(res$avg_log2FC >= 0.378 & res$p_val_adj <= 0.05, '#b2182b', 
                  ifelse(res$avg_log2FC <= -0.378 & res$p_val_adj <= 0.05, '#2166ac',
                         'grey70'))
  names(keyval)[keyval == '#b2182b'] <- 'high'
  names(keyval)[keyval == '#2166ac'] <- 'mid'
  names(keyval)[keyval == 'grey70'] <- 'low'
  #my_xmin = min(res$avg_log2FC, na.rm = T)
  #my_xmax = max(res$avg_log2FC, na.rm = T)
  my_ymax = max(-log10(res$p_val_adj), na.rm = T)
  p1 = EnhancedVolcano(res, lab = lab, 
                       x = 'avg_log2FC', y = 'p_val_adj', 
                       xlab = bquote(~Log[2]~ 'fold change'),
                       ylab = bquote(~-Log[10]~ 'adjusted P'),
                       pCutoff = 0.05, FCcutoff = 0.378, pointSize = 1.5,
                       xlim = xlim,
                       ylim = c(0, my_ymax*1.1),
                       labSize = 4.0, colAlpha = 0.6, legendPosition = 'none',
                       subtitle = '',
                       title = title,
                       caption = '',
                       legendLabSize = 12, legendIconSize = 4.0,
                       titleLabSize = 24,
                       drawConnectors = T,
                       colCustom = keyval,
                       gridlines.major = F,
                       gridlines.minor = F,
                       max.overlaps = 10)
  up_gene_count = paste(dim(up_genes)[1])
  down_gene_count = paste(dim(down_genes)[1])
  p1 = p1 + annotate('text', label = up_gene_count, x = my_xmax, y = my_ymax*1.1, size = 7, 
                     color = '#b2182b', fontface = 'bold') + 
    annotate('text', label = down_gene_count, x = my_xmin, y = my_ymax*1.1, size = 7, 
             color = '#2166ac', fontface = 'bold') +
    theme(plot.margin = margin(t = 0, r = 10, b = -20, l = 10),
          plot.title = element_text(hjust = 0.5, vjust = 0))
  return(p1)
}

gene_voilin_hip_rv <- function(gene, seurat_obj){
  countdata <- FetchData(object = seurat_obj, slot = "data", vars = gene) %>% rownames_to_column("Cell")
  metadata <- seurat_obj@meta.data %>% 
    select(Group, label.main_new) %>% rownames_to_column("Cell")
  metadata <- merge(metadata, countdata, by = "Cell")
  colnames(metadata)[4] = "Expression"
  metadata = metadata[metadata$Expression > 0,]
  
  #only keep cells with >50 cells
  keep <- as.data.frame(table(metadata$label.main_new))
  keep <- as.character(keep[keep$Freq >= 50, ]$Var1)
  
  metadata <- metadata[metadata$label.main_new %in% keep,]
  if(nrow(metadata) == 0){stop("No cells have enough expression")}
  metadata <- metadata[metadata$Group %in% c("Rv", "Hip1"),]
  metadata$Group <- factor(metadata$Group, levels = c("Rv", "Hip1"))
  p1 <- ggplot(metadata, aes(x = Group, y = Expression, fill = Group)) + 
    geom_violin(fill = "white") +
    geom_jitter(size = 0.01, alpha = 0.8) + 
    labs(title = gene) +
    facet_wrap(.~label.main_new, scales = "free") +
    stat_summary(fun = "mean", geom='point', size = 10, colour = "red", shape = 95) +
    stat_compare_means(label = "p.signif", vjust = 1, label.x = 1.4, size = 7, hide.ns = T) +
    theme_Publication() + 
    theme(legend.position = "None", 
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 14))
  return(p1)
}

gene_voilin_hip_rv_macro <- function(gene, seurat_obj){
  countdata <- FetchData(object = seurat_obj, slot = "data", vars = gene) %>% rownames_to_column("Cell")
  metadata <- seurat_obj@meta.data %>% 
    select(Group, label.subset) %>% rownames_to_column("Cell")
  metadata <- merge(metadata, countdata, by = "Cell")
  colnames(metadata)[4] = "Expression"
  metadata = metadata[metadata$Expression > 0,]
  
  #only keep cells with >50 cells
  keep <- as.data.frame(table(metadata$label.subset))
  keep <- as.character(keep[keep$Freq >= 50, ]$Var1)
  
  metadata <- metadata[metadata$label.subset %in% keep,]
  if(nrow(metadata) == 0){stop("No cells have enough expression")}
  metadata <- metadata[metadata$Group %in% c("Rv", "Hip1"),]
  metadata$Group <- factor(metadata$Group, levels = c("Rv", "Hip1"))
  p1 <- ggplot(metadata, aes(x = Group, y = Expression, fill = Group)) + 
    geom_violin(fill = "white") +
    geom_jitter(size = 0.01, alpha = 0.8) + 
    labs(title = gene) +
    facet_wrap(.~label.subset, scales = "free") +
    stat_summary(fun = "mean", geom='point', size = 10, colour = "red", shape = 95) +
    stat_compare_means(label = "p.signif", vjust = 1, label.x = 1.4, size = 7, hide.ns = T) +
    theme_Publication() + 
    theme(legend.position = "None", 
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 14))
  return(p1)
}

gene_voilin_3grps <- function(gene, seurat_obj){
  countdata <- FetchData(object = seurat_obj, slot = "data", vars = gene) %>% rownames_to_column("Cell")
  metadata <- seurat_obj@meta.data %>% 
    select(Group, label.main_new) %>% rownames_to_column("Cell")
  metadata <- merge(metadata, countdata, by = "Cell")
  colnames(metadata)[4] = "Expression"
  metadata = metadata[metadata$Expression > 0,]
  comparison <- list(c("Rv", "Hip1"))
  metadata <- metadata[metadata$Group != "Uninfected",]
  p1 <- ggplot(metadata, aes(x = Group, y = Expression, fill = Group)) + 
    geom_violin() +
    labs(title = gene) +
    facet_wrap(.~label.main_new, nrow = 1) +
    stat_summary(fun = "mean", geom='point', size = 2, colour = "black", show.legend = F) +
    #scale_fill_manual(values = c("#4daf4a", "#377eb8", "#e41a1c"), drop = F) +
    scale_fill_manual(values = c("#377eb8", "#e41a1c")) +
    stat_compare_means(label = "p.signif", vjust = 0.5, hide.ns = T, comparisons = comparison, tip.length = 0.01) +
    theme_Publication() + 
    theme(legend.position = "None", 
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 14))
  return(p1)
}

gene_voilin_2grps_split <- function(gene, seurat_obj){
  countdata <- FetchData(object = seurat_obj, slot = "data", vars = gene) %>% rownames_to_column("Cell")
  metadata <- seurat_obj@meta.data %>% 
    select(Group, label.main_new) %>% rownames_to_column("Cell")
  metadata <- merge(metadata, countdata, by = "Cell")
  colnames(metadata)[4] = "Expression"
  metadata = metadata[metadata$Expression > 0,]
  comparison <- list(c("Rv", "Hip1"))
  metadata <- metadata[metadata$Group != "Uninfected",]
  ymax <- max(metadata$Expression) + 0.3
  plotlist <- lapply(unique(metadata$label.main_new), function(x){
    p1 <- ggplot(data = metadata[metadata$label.main_new == x,], aes(x = Group, y = Expression, fill = Group)) + 
      geom_violin() +
      labs(title = gene) +
      facet_wrap(.~label.main_new, nrow = 1) +
      stat_summary(fun = "mean", geom='point', size = 2, colour = "black", show.legend = F) +
      #scale_fill_manual(values = c("#4daf4a", "#377eb8", "#e41a1c"), drop = F) +
      scale_fill_manual(values = c("#377eb8", "#e41a1c")) +
      scale_y_continuous(limits = c(min(metadata$Expression), ymax)) +
      stat_compare_means(label = "p.signif", vjust = 0.5, hide.ns = T, comparisons = comparison, tip.length = 0.01) +
      theme_Publication() + 
      theme(legend.position = "None", 
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            axis.title.x = element_blank(),
            plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.text = element_text(size = 14))
  })
  names(plotlist) <- unique(metadata$label.main_new)
  return(plotlist)
}

gene_voilin_subset_2grps_split <- function(gene, seurat_obj){
  countdata <- FetchData(object = seurat_obj, slot = "data", vars = gene) %>% rownames_to_column("Cell")
  metadata <- seurat_obj@meta.data %>% 
    select(Group, label.subset) %>% rownames_to_column("Cell")
  metadata <- merge(metadata, countdata, by = "Cell")
  colnames(metadata)[4] = "Expression"
  metadata = metadata[metadata$Expression > 0,]
  comparison <- list(c("Rv", "Hip1"))
  metadata <- metadata[metadata$Group != "Uninfected",]
  ymax <- max(metadata$Expression) + 0.3
  plotlist <- lapply(unique(metadata$label.subset), function(x){
    p1 <- ggplot(data = metadata[metadata$label.subset == x,], aes(x = Group, y = Expression, fill = Group)) + 
      geom_violin() +
      labs(title = gene) +
      facet_wrap(.~label.subset, nrow = 1) +
      stat_summary(fun = "mean", geom='point', size = 2, colour = "black", show.legend = F) +
      #scale_fill_manual(values = c("#4daf4a", "#377eb8", "#e41a1c"), drop = F) +
      scale_fill_manual(values = c("#377eb8", "#e41a1c")) +
      scale_y_continuous(limits = c(min(metadata$Expression), ymax)) +
      stat_compare_means(label = "p.signif", vjust = 0.5, hide.ns = T, comparisons = comparison, tip.length = 0.01) +
      theme_Publication() + 
      theme(legend.position = "None", 
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            axis.title.x = element_blank(),
            plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.text = element_text(size = 14))
  })
  names(plotlist) <- unique(metadata$label.subset)
  return(plotlist)
}

ref.genelist = function(genes){
  genenames = genes$GeneID
  genenames <- stringr::str_to_upper(genenames)
  #annotate gene names to entrez_id
  ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
  
  entrez_names = getBM(attributes = c('external_gene_name', 'entrezgene_id'), 
                       filters = 'external_gene_name',
                       values = genenames, mart = ensembl)
  
  genelist = data.frame('external_gene_name' = genes$GeneID, 
                        'foldchange' = genes$avg_log2FC)
  
  genelist = merge(genelist, entrez_names, by = 'external_gene_name')
  genelist = genelist[!is.na(genelist$entrezgene_id),]
  kegg_genes = genelist[,2]
  names(kegg_genes) = as.character(genelist[,3])
  return(kegg_genes)
}

pathview_meta <- function (gene.data = NULL, cpd.data = NULL, pathway.id, species = "hsa", 
                           kegg.dir = ".", cpd.idtype = "kegg", gene.idtype = "entrez", 
                           gene.annotpkg = NULL, min.nnodes = 3, kegg.native = TRUE, 
                           map.null = TRUE, expand.node = FALSE, split.group = FALSE, 
                           map.symbol = TRUE, map.cpdname = TRUE, node.sum = "sum", 
                           discrete = list(gene = FALSE, cpd = FALSE), limit = list(gene = 1, 
                                                                                    cpd = 1), bins = list(gene = 10, cpd = 10), both.dirs = list(gene = T, 
                                                                                                                                                 cpd = T), trans.fun = list(gene = NULL, cpd = NULL), 
                           low = list(gene = "green", cpd = "blue"), mid = list(gene = "gray", 
                                                                                cpd = "gray"), high = list(gene = "red", cpd = "yellow"), 
                           na.col = "transparent", ...) 
{
  dtypes = !is.null(gene.data) + !is.null(cpd.data)
  cond0 = dtypes == 1 & is.numeric(limit) & length(limit) > 
    1
  if (cond0) {
    if (limit[1] != limit[2] & is.null(names(limit))) 
      limit = list(gene = limit[1:2], cpd = limit[1:2])
  }
  if (is.null(trans.fun)) 
    trans.fun = list(gene = NULL, cpd = NULL)
  arg.len2 = c("discrete", "limit", "bins", "both.dirs", "trans.fun", 
               "low", "mid", "high")
  for (arg in arg.len2) {
    obj1 = eval(as.name(arg))
    if (length(obj1) == 1) 
      obj1 = rep(obj1, 2)
    if (length(obj1) > 2) 
      obj1 = obj1[1:2]
    obj1 = as.list(obj1)
    ns = names(obj1)
    if (length(ns) == 0 | !all(c("gene", "cpd") %in% ns)) 
      names(obj1) = c("gene", "cpd")
    assign(arg, obj1)
  }
  if (is.character(gene.data)) {
    gd.names = gene.data
    gene.data = rep(1, length(gene.data))
    names(gene.data) = gd.names
    both.dirs$gene = FALSE
    ng = length(gene.data)
    nsamp.g = 1
  }
  else if (!is.null(gene.data)) {
    if (length(dim(gene.data)) == 2) {
      gd.names = rownames(gene.data)
      ng = nrow(gene.data)
      nsamp.g = 2
    }
    else if (is.numeric(gene.data) & is.null(dim(gene.data))) {
      gd.names = names(gene.data)
      ng = length(gene.data)
      nsamp.g = 1
    }
    else stop("wrong gene.data format!")
  }
  else if (is.null(cpd.data)) {
    stop("gene.data and cpd.data are both NULL!")
  }
  gene.idtype = toupper(gene.idtype)
  data(bods)
  if (species != "ko") {
    species.data = kegg.species.code(species, na.rm = T, 
                                     code.only = FALSE)
  }
  else {
    species.data = c(kegg.code = "ko", entrez.gnodes = "0", 
                     kegg.geneid = "K01488", ncbi.geneid = NA, ncbi.proteinid = NA, 
                     uniprot = NA)
    gene.idtype = "KEGG"
    msg.fmt = "Only KEGG ortholog gene ID is supported, make sure it looks like \"%s\"!"
    msg = sprintf(msg.fmt, species.data["kegg.geneid"])
    message("Note: ", msg)
  }
  if (length(dim(species.data)) == 2) {
    message("Note: ", "More than two valide species!")
    species.data = species.data[1, ]
  }
  species = species.data["kegg.code"]
  entrez.gnodes = species.data["entrez.gnodes"] == 1
  if (is.na(species.data["ncbi.geneid"])) {
    if (!is.na(species.data["kegg.geneid"])) {
      msg.fmt = "Mapping via KEGG gene ID (not Entrez) is supported for this species,\nit looks like \"%s\"!"
      msg = sprintf(msg.fmt, species.data["kegg.geneid"])
      message("Note: ", msg)
    }
    else {
      stop("This species is not annotated in KEGG!")
    }
  }
  if (is.null(gene.annotpkg)) 
    gene.annotpkg = bods[match(species, bods[, 3]), 1]
  if (length(grep("ENTREZ|KEGG|NCBIPROT|UNIPROT", gene.idtype)) < 
      1 & !is.null(gene.data)) {
    if (is.na(gene.annotpkg)) 
      stop("No proper gene annotation package available!")
    if (!gene.idtype %in% gene.idtype.bods[[species]]) 
      stop("Wrong input gene ID type!")
    gene.idmap = id2eg(gd.names, category = gene.idtype, 
                       pkg.name = gene.annotpkg, unique.map = F)
    gene.data = mol.sum(gene.data, gene.idmap)
    gene.idtype = "ENTREZ"
  }
  if (gene.idtype != "KEGG" & !entrez.gnodes & !is.null(gene.data)) {
    id.type = gene.idtype
    if (id.type == "ENTREZ") 
      id.type = "ENTREZID"
    kid.map = names(species.data)[-c(1:2)]
    kid.types = names(kid.map) = c("KEGG", "ENTREZID", "NCBIPROT", 
                                   "UNIPROT")
    kid.map2 = gsub("[.]", "-", kid.map)
    kid.map2["UNIPROT"] = "up"
    if (is.na(kid.map[id.type])) 
      stop("Wrong input gene ID type for the species!")
    message("Info: Getting gene ID data from KEGG...")
    gene.idmap = keggConv(kid.map2[id.type], species)
    message("Info: Done with data retrieval!")
    kegg.ids = gsub(paste(species, ":", sep = ""), "", names(gene.idmap))
    in.ids = gsub(paste0(kid.map2[id.type], ":"), "", gene.idmap)
    gene.idmap = cbind(in.ids, kegg.ids)
    gene.data = mol.sum(gene.data, gene.idmap)
    gene.idtype = "KEGG"
  }
  if (is.character(cpd.data)) {
    cpdd.names = cpd.data
    cpd.data = rep(1, length(cpd.data))
    names(cpd.data) = cpdd.names
    both.dirs$cpd = FALSE
    ncpd = length(cpd.data)
  }
  else if (!is.null(cpd.data)) {
    if (length(dim(cpd.data)) == 2) {
      cpdd.names = rownames(cpd.data)
      ncpd = nrow(cpd.data)
    }
    else if (is.numeric(cpd.data) & is.null(dim(cpd.data))) {
      cpdd.names = names(cpd.data)
      ncpd = length(cpd.data)
    }
    else stop("wrong cpd.data format!")
  }
  if (length(grep("kegg", cpd.idtype)) < 1 & !is.null(cpd.data)) {
    data(rn.list)
    cpd.types = c(names(rn.list), "name")
    cpd.types = tolower(cpd.types)
    cpd.types = cpd.types[-grep("kegg", cpd.types)]
    if (!tolower(cpd.idtype) %in% cpd.types) 
      stop("Wrong input cpd ID type!")
    cpd.idmap = cpd2kegg(cpdd.names, in.type = cpd.idtype)
    cpd.data = mol.sum(cpd.data, cpd.idmap)
  }
  warn.fmt = "Parsing %s file failed, please check the file!"
  if (length(grep(species, pathway.id)) > 0) {
    pathway.name = pathway.id
    pathway.id = gsub(species, "", pathway.id)
  }
  else pathway.name = paste(species, pathway.id, sep = "")
  kfiles = list.files(path = kegg.dir, pattern = "[.]xml|[.]png")
  npath = length(pathway.id)
  out.list = list()
  tfiles.xml = paste(pathway.name, "xml", sep = ".")
  tfiles.png = paste(pathway.name, "png", sep = ".")
  if (kegg.native) 
    ttype = c("xml", "png")
  else ttype = "xml"
  xml.file <- paste(kegg.dir, "/", tfiles.xml, sep = "")
  for (i in 1:npath) {
    if (kegg.native) 
      tfiles = c(tfiles.xml[i], tfiles.png[i])
    else tfiles = tfiles.xml[i]
    if (!all(tfiles %in% kfiles)) {
      dstatus = download.kegg(pathway.id = pathway.id[i], 
                              species = species, kegg.dir = kegg.dir, file.type = ttype)
      if (dstatus == "failed") {
        warn.fmt = "Failed to download KEGG xml/png files, %s skipped!"
        warn.msg = sprintf(warn.fmt, pathway.name[i])
        message("Warning: ", warn.msg)
        return(invisible(0))
      }
    }
    if (kegg.native) {
      node.data = try(node.info(xml.file[i]), silent = T)
      if (class(node.data)[1] == "try-error") {
        warn.msg = sprintf(warn.fmt, xml.file[i])
        message("Warning: ", warn.msg)
        return(invisible(0))
      }
      node.type = c("gene", "enzyme", "compound", "ortholog")
      sel.idx = node.data$type %in% node.type
      nna.idx = !is.na(node.data$x + node.data$y + node.data$width + 
                         node.data$height)
      sel.idx = sel.idx & nna.idx
      if (sum(sel.idx) < min.nnodes) {
        warn.fmt = "Number of mappable nodes is below %d, %s skipped!"
        warn.msg = sprintf(warn.fmt, min.nnodes, pathway.name[i])
        message("Warning: ", warn.msg)
        return(invisible(0))
      }
      node.data = lapply(node.data, "[", sel.idx)
    }
    else {
      gR1 = try(parseKGML2Graph2(xml.file[i], genes = F, 
                                 expand = expand.node, split.group = split.group), 
                silent = T)
      node.data = try(node.info(gR1), silent = T)
      if (class(node.data)[1] == "try-error") {
        warn.msg = sprintf(warn.fmt, xml.file[i])
        message("Warning: ", warn.msg)
        return(invisible(0))
      }
    }
    if (species == "ko") 
      gene.node.type = "ortholog"
    else gene.node.type = "gene"
    if ((!is.null(gene.data) | map.null) & sum(node.data$type == 
                                               gene.node.type) > 1) {
      plot.data.gene = node.map(gene.data, node.data, 
                                node.types = gene.node.type, node.sum = node.sum, 
                                entrez.gnodes = entrez.gnodes)
      kng = plot.data.gene$kegg.names
      kng.char = gsub("[0-9]", "", unlist(kng))
      if (any(kng.char > "")) 
        entrez.gnodes = FALSE
      if (map.symbol & species != "ko" & entrez.gnodes) {
        if (is.na(gene.annotpkg)) {
          warn.fmt = "No annotation package for the species %s, gene symbols not mapped!"
          warn.msg = sprintf(warn.fmt, species)
          message("Warning: ", warn.msg)
        }
        else {
          plot.data.gene$labels = eg2id(as.character(plot.data.gene$kegg.names), 
                                        category = "SYMBOL", pkg.name = gene.annotpkg)[, 
                                                                                       2]
          mapped.gnodes = rownames(plot.data.gene)
          node.data$labels[mapped.gnodes] = plot.data.gene$labels
        }
      }
      cols.ts.gene = node.color(plot.data.gene, limit$gene, 
                                bins$gene, both.dirs = both.dirs$gene, trans.fun = trans.fun$gene, 
                                discrete = discrete$gene, low = low$gene, mid = mid$gene, 
                                high = high$gene, na.col = na.col)
    }
    else plot.data.gene = cols.ts.gene = NULL
    if ((!is.null(cpd.data) | map.null) & sum(node.data$type == 
                                              "compound") > 1) {
      plot.data.cpd = node.map(cpd.data, node.data, node.types = "compound", 
                               node.sum = node.sum)
      if (map.cpdname & !kegg.native) {
        plot.data.cpd$labels = cpdkegg2name(plot.data.cpd$labels)[, 
                                                                  2]
        mapped.cnodes = rownames(plot.data.cpd)
        node.data$labels[mapped.cnodes] = plot.data.cpd$labels
      }
      cols.ts.cpd = node.color(plot.data.cpd, limit$cpd, 
                               bins$cpd, both.dirs = both.dirs$cpd, trans.fun = trans.fun$cpd, 
                               discrete = discrete$cpd, low = low$cpd, mid = mid$cpd, 
                               high = high$cpd, na.col = na.col)
    }
    else plot.data.cpd = cols.ts.cpd = NULL
    if (kegg.native) {
      pv.pars = keggview.native(plot.data.gene = plot.data.gene, 
                                cols.ts.gene = cols.ts.gene, plot.data.cpd = plot.data.cpd, 
                                cols.ts.cpd = cols.ts.cpd, node.data = node.data, 
                                pathway.name = pathway.name[i], kegg.dir = kegg.dir, 
                                limit = limit, bins = bins, both.dirs = both.dirs, 
                                discrete = discrete, low = low, mid = mid, high = high, 
                                na.col = na.col, ...)
    }
    else {
      pv.pars = keggview.graph(plot.data.gene = plot.data.gene, 
                               cols.ts.gene = cols.ts.gene, plot.data.cpd = plot.data.cpd, 
                               cols.ts.cpd = cols.ts.cpd, node.data = node.data, 
                               path.graph = gR1, pathway.name = pathway.name[i], 
                               map.cpdname = map.cpdname, split.group = split.group, 
                               limit = limit, bins = bins, both.dirs = both.dirs, 
                               discrete = discrete, low = low, mid = mid, high = high, 
                               na.col = na.col, ...)
    }
    plot.data.gene = cbind(plot.data.gene, cols.ts.gene)
    if (!is.null(plot.data.gene)) {
      cnames = colnames(plot.data.gene)[-(1:8)]
      nsamp = length(cnames)/2
      if (nsamp > 1) {
        cnames[(nsamp + 1):(2 * nsamp)] = paste(cnames[(nsamp + 
                                                          1):(2 * nsamp)], "col", sep = ".")
      }
      else cnames[2] = "mol.col"
      colnames(plot.data.gene)[-(1:8)] = cnames
    }
    plot.data.cpd = cbind(plot.data.cpd, cols.ts.cpd)
    if (!is.null(plot.data.cpd)) {
      cnames = colnames(plot.data.cpd)[-(1:8)]
      nsamp = length(cnames)/2
      if (nsamp > 1) {
        cnames[(nsamp + 1):(2 * nsamp)] = paste(cnames[(nsamp + 
                                                          1):(2 * nsamp)], "col", sep = ".")
      }
      else cnames[2] = "mol.col"
      colnames(plot.data.cpd)[-(1:8)] = cnames
    }
    out.list[[i]] = list(plot.data.gene = plot.data.gene, 
                         plot.data.cpd = plot.data.cpd)
  }
  if (npath == 1) 
    out.list = out.list[[1]]
  else names(out.list) = pathway.name
  return(plot.data.gene)
}


###################
ggmaplot_cus <- function (data, fdr = 0.05, fc = 1.5, genenames = NULL, detection_call = NULL, size = NULL, alpha = 1, seed = 42, font.label = c(12, "plain", "black"), label.rectangle = FALSE, palette = c("#B31B21", "#1465AC", "darkgray"), top = 15, select.top.method = c("padj",  "fc"), label.select = NULL, main = NULL, xlab = "Log2 mean expression",  ylab = "Log2 fold change", ggtheme = theme_classic(), ...) 
{
  if (!base::inherits(data, c("matrix", "data.frame", "DataFrame", 
                              "DE_Results", "DESeqResults"))) 
    stop("data must be an object of class matrix, data.frame, DataFrame, DE_Results or DESeqResults")
  if (!is.null(detection_call)) {
    if (nrow(data) != length(detection_call)) 
      stop("detection_call must be a numeric vector of length = nrow(data)")
  }
  else if ("detection_call" %in% colnames(data)) {
    detection_call <- as.vector(data$detection_call)
  }
  else detection_call = rep(1, nrow(data))
  if (is.null(list(...)$legend)) 
    legend <- c(0.12, 0.9)
  is.basemean.logged <- "baseMeanLog2" %in% colnames(data)
  if (is.basemean.logged) {
    data$baseMean <- data$baseMeanLog2
  }
  else if ("baseMean" %in% colnames(data)) {
    data$baseMean <- log2(data$baseMean + 1)
  }
  ss <- base::setdiff(c("baseMean", "log2FoldChange", "padj"), 
                      colnames(data))
  if (length(ss) > 0) 
    stop("The colnames of data must contain: ", paste(ss, 
                                                      collapse = ", "))
  if (is.null(genenames)) 
    genenames <- rownames(data)
  else if (length(genenames) != nrow(data)) 
    stop("genenames should be of length nrow(data).")
  sig <- rep(3, nrow(data))
  sig[which(data$padj <= fdr & data$log2FoldChange < 0 & abs(data$log2FoldChange) >= 
              log2(fc) & detection_call == 1)] = 2
  sig[which(data$padj <= fdr & data$log2FoldChange > 0 & abs(data$log2FoldChange) >= 
              log2(fc) & detection_call == 1)] = 1
  data <- data.frame(name = genenames, mean = data$baseMean, 
                     lfc = data$log2FoldChange, padj = data$padj, sig = sig)
  . <- NULL
  data$sig <- as.factor(data$sig)
  .lev <- .levels(data$sig) %>% as.numeric()
  palette <- palette[.lev]
  new.levels <- c(paste0("Up: ", sum(sig == 1)), paste0("Down: ", 
                                                        sum(sig == 2)), "NS") %>% .[.lev]
  data$sig <- factor(data$sig, labels = new.levels)
  select.top.method <- match.arg(select.top.method)
  if (select.top.method == "padj") 
    data <- data[order(data$padj), ]
  else if (select.top.method == "fc") 
    data <- data[order(abs(data$lfc), decreasing = TRUE), 
    ]
  complete_data <- stats::na.omit(data)
  labs_data <- subset(complete_data, padj <= fdr & name != 
                        "" & abs(lfc) >= log2(fc))
  labs_data <- utils::head(labs_data, top)
  if (!is.null(label.select)) {
    selected_labels <- complete_data %>% subset(complete_data$name %in% 
                                                  label.select, drop = FALSE)
    labs_data <- dplyr::bind_rows(labs_data, selected_labels) %>% 
      dplyr::distinct(.data$name, .keep_all = TRUE)
  }
  font.label <- .parse_font(font.label)
  font.label$size <- ifelse(is.null(font.label$size), 12, 
                            font.label$size)
  font.label$color <- ifelse(is.null(font.label$color), "black", 
                             font.label$color)
  font.label$face <- ifelse(is.null(font.label$face), "plain", 
                            font.label$face)
  mean <- lfc <- sig <- name <- padj <- NULL
  p <- ggplot(data, aes(x = mean, y = lfc)) + geom_point(aes(color = sig), 
                                                         size = size, alpha = alpha)
  max.overlaps = getOption("ggrepel.max.overlaps", default = Inf)
  if (label.rectangle) {
    p <- p + ggrepel::geom_label_repel(data = labs_data, 
                                       mapping = aes(label = name), box.padding = unit(0.35, 
                                                                                       "lines"), point.padding = unit(0.3, "lines"), 
                                       force = 1, seed = seed, fontface = font.label$face, 
                                       size = font.label$size/3, color = font.label$color, 
                                       max.overlaps = max.overlaps)
  }
  else {
    p <- p + ggrepel::geom_text_repel(data = labs_data, 
                                      mapping = aes(label = name), box.padding = unit(0.35, 
                                                                                      "lines"), point.padding = unit(0.3, "lines"), 
                                      force = 1, seed = seed, fontface = font.label$face, 
                                      size = font.label$size/3, color = font.label$color, 
                                      max.overlaps = max.overlaps)
  }
  p <- p + scale_x_continuous(breaks = seq(0, max(data$mean), 
                                           2)) + labs(x = xlab, y = ylab, title = main, color = "") + 
    geom_hline(yintercept = c(0, -log2(fc), log2(fc)), linetype = c(1, 
                                                                    2, 2), color = c("black", "black", "black"))
  p <- ggpar(p, palette = palette, ggtheme = ggtheme, ...)
  p
}

genesigs.umap.dens.Aucell <- function(mygenes, pathwayname, seuratobj = lung9s3){
  library(AUCell)
  library(GSEABase)
  #Aucell scores
  exprMatrix <- seuratobj[["RNA"]]$counts
  geneSets <- list(mygenes)
  names(geneSets) <- pathwayname
  cells_AUC <- AUCell_run(exprMatrix, geneSets, BPPARAM=BiocParallel::MulticoreParam(8))
  scores <- getAUC(cells_AUC)[rownames(cells_AUC),]
  pathtitle <- paste(pathwayname, "_AuCell", sep = "")
  scores <- data.frame(scores)
  colnames(scores) <- pathtitle
  seuratobj <- AddMetaData(seuratobj, metadata = scores)
  
  p1 <- FeaturePlot(subset(seuratobj, subset = Group != "Uninfected"), features = pathtitle, split.by = "Group", pt.size = 0.5, alpha = 0.8) + plot_layout(guides = "collect") & scale_color_gradient(low = "blue", high = "red") & theme(legend.position = "bottom")
  
  #Density plots
  meta <- seuratobj@meta.data
  meta <- split(meta, f = meta$label.main_new)
  
  out <- lapply(meta, function(x){
    df <- x[c(grep("Group", colnames(x), ignore.case = T), grep("AuCell", colnames(x), ignore.case = T))]
    colnames(df) <- c("Group", "Pathway")
    pval <- t.test(df[df$Group == "hip1 Mut",2], df[df$Group == "WT Mtb",2], alternative = "greater")$p.value
    pval <- format.pval(pv = pval, digits = 2, eps = 0.0001, nsmall = 3)
    xmax <- max(density(df[df$Group != "Uninfected",]$Pathway)$x)
    ymax <- max(density(df[df$Group != "Uninfected",]$Pathway)$y)
    ggplot(df[df$Group != "Uninfected",], aes(x = Pathway, fill = Group)) + geom_density(alpha = 0.5) + labs(title = unique(x$label.main_new)) + annotate("text", x = xmax/2, y = ymax, label = paste("p-value: ", pval,sep = "")) + theme_Publication() + scale_fill_Publication() + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5, vjust = -1), plot.margin = margin(0,0,0,0)) 
  })
  p2 <- wrap_plots(out) + plot_layout(nrow = 4, ncol = 4, guides = "collect") & theme(legend.position = "bottom")
  
  out <- lapply(meta, function(x){
    df <- x[c(grep("Group", colnames(x), ignore.case = T), grep("AuCell", colnames(x), ignore.case = T))]
    colnames(df) <- c("Group", "Pathway")
    ggplot(df[df$Group != "Uninfected",], aes(x = Group, y = Pathway, fill = Group)) + geom_violin(alpha = 0.5) + labs(title = unique(x$label.main_new)) + stat_compare_means(label = "p.signif", label.x = 1.465, size = 6, vjust = 1) + theme_Publication() + scale_fill_Publication() + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5, vjust = -1), plot.margin = margin(0,0,0,0)) 
  })
  p3 <- wrap_plots(out) + plot_layout(nrow = 3, ncol = 5, guides = "collect") & theme(legend.position = "bottom") 
  out <- list(p1, p2, p3)
}

genesigs.umap.dens.Sipsic <- function(mygenes, pathwayname, seuratobj = lung9s3){
  library(SiPSiC)
  #Sipsic scores
  seuratobj_sce <- as.SingleCellExperiment(seuratobj)
  scoresAndIndices <- getPathwayScores(counts(seuratobj_sce), mygenes)
  pathwayScoresOfCells <- scoresAndIndices$pathwayScores
  pathtitle <- paste(pathwayname, "_SipSic", sep = "")
  pathwayScoresOfCells <- data.frame(pathwayScoresOfCells)
  colnames(pathwayScoresOfCells) <- pathtitle
  seuratobj <- AddMetaData(seuratobj, metadata = pathwayScoresOfCells)
  
  p1 <- FeaturePlot(subset(seuratobj, subset = Group != "Uninfected"), features = pathtitle, split.by = "Group", pt.size = 0.5, alpha = 0.8) + plot_layout(guides = "collect") & scale_color_gradient(low = "black", high = "yellow", limits = c(0,0.5)) & theme(legend.position = "bottom")
  
  #Density plots
  meta <- seuratobj@meta.data
  meta <- split(meta, f = meta$label.main_new)
  
  out <- lapply(meta, function(x){
    df <- x[c(grep("Group", colnames(x), ignore.case = T), grep("SipSic", colnames(x), ignore.case = T))]
    colnames(df) <- c("Group", "Pathway")
    pval <- t.test(df[df$Group == "hip1 Mut",2], df[df$Group == "WT Mtb",2], alternative = "greater")$p.value
    pval <- format.pval(pv = pval, digits = 2, eps = 0.0001, nsmall = 3)
    xmax <- max(density(df[df$Group != "Uninfected",]$Pathway)$x)
    ymax <- max(density(df[df$Group != "Uninfected",]$Pathway)$y)
    ggplot(df[df$Group != "Uninfected",], aes(x = Pathway, fill = Group)) + geom_density(alpha = 0.5) + labs(title = unique(x$label.main_new)) + annotate("text", x = xmax/2, y = ymax, label = paste("p-value: ", pval,sep = "")) + theme_Publication() + scale_fill_Publication() + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5, vjust = -1), plot.margin = margin(0,0,0,0)) 
  })
  p2 <- wrap_plots(out) + plot_layout(nrow = 4, ncol = 4, guides = "collect") & theme(legend.position = "bottom")
  
  out <- lapply(meta, function(x){
    df <- x[c(grep("Group", colnames(x), ignore.case = T), grep("SipSic", colnames(x), ignore.case = T))]
    colnames(df) <- c("Group", "Pathway")
    ggplot(df[df$Group != "Uninfected",], aes(x = Group, y = Pathway, fill = Group)) + geom_violin(alpha = 0.5) + labs(title = unique(x$label.main_new)) + stat_compare_means(label = "p.signif", label.x = 1.465, size = 6, vjust = 1) + theme_Publication() + scale_fill_Publication() + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5, vjust = -1), plot.margin = margin(0,0,0,0)) 
  })
  p3 <- wrap_plots(out) + plot_layout(nrow = 3, ncol = 5, guides = "collect") & theme(legend.position = "bottom") 
  out <- list(p1, p2, p3)
}

maplot_data <- function(seurat_obj, subsetmain, toptable_dat){
  test <- AverageExpression(subset(seurat_obj, subset = label.main_new == subsetmain), group.by = "orig.ident")[[1]]
  test <- as.data.frame(test)
  test$log2exp <- log2(test$all + 1)
  test <- test %>% rownames_to_column("GeneID")
  plotdat_merge <- merge(toptable_dat, test, by = "GeneID")
  plotdat_merge <- plotdat_merge[c("GeneID", "avg_log2FC", "p_val_adj", "log2exp")]
  #topTables_global <- topTables_global %>% column_to_rownames("GeneID")
  colnames(plotdat_merge) <- c("GeneID","log2FoldChange", "padj", "baseMeanLog2")
  
  genes_up <- toptable_dat[toptable_dat$avg_log2FC >= 0.378 & toptable_dat$p_val_adj <= 0.05,]
  genes_down <- toptable_dat[toptable_dat$avg_log2FC <= -0.378 & toptable_dat$p_val_adj < 0.05,]
  
  set.seed(101010)
  data <- ggmaplot(plotdat_merge, fdr = 0.05, fc = 1.3, size = 2.5, palette = c("red", "blue", "darkgray"), genenames = as.vector(plotdat_merge$GeneID), label.select = c(rownames(genes_up), rownames(genes_down)), legend = "top", top = 20, font.label = c("bold", 12), label.rectangle = F, font.legend = "bold", font.main = "bold", ggtheme = ggplot2::theme_minimal()) + scale_y_continuous(limits = c(-4, 4)) + theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12), axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), legend.text = element_text(size = 16)) +  guides(colour = guide_legend(override.aes = list(size=10)))
  data <- data$data
  data$size <- ifelse(data$sig == "NS", 2, 3.3)
  data$alpha <- ifelse(data$sig == "NS", 0.99, 1)
  data$lfc[data$lfc >= 3 | data$lfc <= -3] <- datawizard::rescale(data$lfc[data$lfc >= 3 | data$lfc <= -3], to = c(-2,2))
  return(data)
}