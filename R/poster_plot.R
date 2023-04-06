# Sys.getenv("GITHUB_PAT")
# Sys.unsetenv("GITHUB_PAT")

library(reticulate)
reticulate::use_condaenv("/opt/anaconda3/envs/basilica-env/bin/python")
devtools::load_all("~/Documents/GitHub/basilica")


library(CNAqc)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(reshape2)




# data

x <- readRDS("~/Documents/GitHub/simbasilica/script_test/simulations/simul.N100.G2.s23.Rds")

dn = x$exp_denovo %>% as.data.frame() %>% rownames
sc = x$exp_fixed %>% as.data.frame() %>% rownames

dn_c = ggsci::pal_nejm()(length(dn))
names(dn_c) = dn

sc_c = ggsci::pal_simpsons()(length(sc))
names(sc_c) = sc

cls = c(sc_c,dn_c)

b = x$exp_exposure
a = rbind(x$exp_denovo %>% as.data.frame(),x$exp_fixed %>% as.data.frame())

my_plot_exposure = function(b,cls){

  ggplot(data = b %>% as.data.frame() %>% mutate(sample = paste0(1:100)) %>%
                      melt() %>% dplyr::rename(Signature = variable),aes(x = sample, y  = value, fill = Signature)) +
  geom_bar(stat = "identity")  + ggplot2::scale_fill_manual(values = cls ) + labs(title = "Expsosure") + theme(

    axis.ticks.x = element_blank(),axis.text.x = element_blank()
  )

}


my_plot_signatures = function(a,cls,levels){

 a %>% as.data.frame() %>% dplyr::mutate(sbs = rownames(a)) %>%
  reshape2::melt() %>%
  as_tibble() %>% dplyr::rename(Var1 = sbs, Var2 = variable) %>%
  mutate(
    substitution = paste0(substr(start = 3, stop = 3, Var2),">",substr(start = 5, stop = 5, Var2)),
    context = paste0(
      substr(start = 1, stop = 1, Var2),
      '_',
      substr(start = 7, stop = 7, Var2)
    )
  )  %>%
  ggplot() +
  geom_bar(aes(value, x = context, fill = Var1), stat = 'identity') +
  facet_grid(factor(Var1, levels= levels) ~ substitution, scales = 'free') +
  my_ggplot_theme() +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank()) +
  scale_fill_manual(values = cls) +
  guides(fill = 'none')  +
  labs(x = '', y = "", title = "Signatures")

}

plot_simulated_data = function(x, cls = c("1"= "forestgreen", "2"= "purple3")){

  x$x[[1]] %>%  dplyr::mutate(group=paste0(x$groups[[1]])) %>%
  reshape2::melt() %>%
  as_tibble() %>%
  mutate(
    substitution = substr(start = 3, stop = 5, variable),
    context = paste0(
      substr(start = 1, stop = 1, variable),
      '_',
      substr(start = 7, stop = 7, variable)
    )
  ) %>%  dplyr::group_by(group) %>%
  dplyr::mutate(n_muts=value/sum(value)) %>%
  ggplot() +
  geom_bar(aes(y = n_muts, x = context, fill = group), stat = 'identity') +
  my_ggplot_theme()  +
  facet_grid(group~ substitution, scales = 'free') + scale_fill_manual(values = cls) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank()
  ) + labs(x = '', y = "Mutation count")  + labs(title = "Data")

}


expos = my_plot_exposure(b,cls)

sig = my_plot_signatures(a,cls,levels = c("SBS1","SBS10a","SBS17b","SBS4","SBS5","SBS6"))

data = plot_simulated_data(x) +
  theme(legend.position = "none")

# plot_inference

x.fit.noreg <- readRDS("~/Documents/GitHub/simbasilica/script_test/simulations/fit_noreg_nogroups_no_private.Rds")


p1 = plot_signatures(x.fit.noreg, Type = 'Catalogue') + labs(title = 'Inferred Catalogue Signatures', x = "", y = "") +
theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.position = "none")

ggsave(p1,filename = "inferred_catalogue_example.png",width = 7,height = 3)

names(cls) = c("SBS1","SBS5","SBS17b","D2","D3","D1")

inf_sign = rbind(as.data.frame(x.fit.noreg$fit$denovo_signatures),as.data.frame(x.fit.noreg$fit$catalogue_signatures)) %>%
  my_plot_signatures(cls,levels = c("SBS1","D1","SBS17b","D2","SBS5","D3")) +
  labs(title = "Inferred Signatures", x= "")

inf_exp = x.fit.noreg$fit$exposure %>% my_plot_exposure(cls) +  labs(title = "Inferred Exposure")

signatures = rbind(x.fit.noreg$fit$catalogue_signatures,x.fit.noreg$fit$denovo_signatures)
exp = x.fit.noreg$fit$exposure

tmb = x$x[[1]] %>% rowSums()

rec_data =  as.matrix(exp*tmb) %*%  as.matrix(signatures) %>%  as.data.frame() %>%
  dplyr::mutate(group=paste0(x$groups[[1]])) %>%
  reshape2::melt() %>%
  as_tibble() %>%
  mutate(
    substitution = substr(start = 3, stop = 5, variable),
    context = paste0(
      substr(start = 1, stop = 1, variable),
      '_',
      substr(start = 7, stop = 7, variable)
    )
  ) %>%  dplyr::group_by(group) %>%
  dplyr::mutate(n_muts=value/sum(value)) %>%
  ggplot() +
  geom_bar(aes(y = n_muts, x = context, fill = group), stat = 'identity') +
  my_ggplot_theme()  +
  facet_grid(group~ substitution, scales = 'free') + scale_fill_manual(values = c("1"= "forestgreen", "2"= "purple3")) +
  theme(
    strip.text.y = element_text(angle = 0),
    axis.text.x = element_text(angle = 90, size = 4)
  ) + labs(x = 'Context', y = "Mutation count", title = "Reconstructed data")


setwd("~/Documents/GitHub/simbasilica")

# scores

cosmic = read.csv("./script_test/COSMIC_v3.3.1_SBS_GRCh38.txt", sep="\t") %>%
  tibble::column_to_rownames(var="Type") %>% t()

input_list = c("SBS1","SBS5","SBS6","SBS7d","SBS33","SBS22","SBS10a","SBS4")

sim_ref = x.fit.noreg %>% plot_similarity_reference(reference = cosmic[input_list,],context = F) +  plot_annotation(title = "De Novo discovery")

reconstr_data =  as.matrix(exp*tmb) %*%  as.matrix(signatures) %>%  as.data.frame()

real_data = x$x[[1]]

cosines = lapply(1:nrow(reconstr_data),function(i){

    tibble(sample = paste0(i), cos =  cosine.vector(reconstr_data[i,],real_data[i,]))

}) %>% bind_rows()

cosine_plot = ggplot(cosines %>% dplyr::mutate(group=paste0(x$groups[[1]])),aes(x = group,y = cos, fill = group)) + geom_boxplot() +
  CNAqc:::my_ggplot_theme() + scale_fill_manual(values = c("1"= "forestgreen", "2"= "purple3")) +
  labs(y = "Cosine Similarity", title = "Data reconstruction")  + theme(
    axis.text.x = element_blank()
  )


ggsave(sig,filename = paste0("input_sign.png"),height = 4,width = 6)
ggsave(expos,filename = paste0("input_expos.png"),height = 3,width = 6)
ggsave(data + plot_annotation("Example of simulated data"),filename = paste0("input_data.png"),height = 4,width = 6)

ggsave(inf_sign,filename = paste0("inf_sign.png"),height = 3,width = 6)
ggsave(inf_exp,filename = paste0("inf_expos.png"),height = 3,width = 6)
ggsave(rec_data,filename = paste0("rec_data.png"),height = 4,width = 6)

ggsave(sim_ref ,filename = paste0("sim_ref.png"),height = 7,width = 8)




















# cnv statistics




plot_data_signatures = function(a,what = "SBS", context = T, cls = NULL){

  a =  a %>% dplyr::mutate(sbs = rownames(a)) %>% as_tibble() %>%
    reshape2::melt()
  library(stringr)

  if(what == 'SBS'){
    a = a %>% dplyr::rename(Var1 = sbs, Var2 = variable) %>%
      mutate(
        substitution = paste0(substr(start = 3, stop = 3, Var2),">",substr(start = 5, stop = 5, Var2)),
        context = paste0(
          substr(start = 1, stop = 1, Var2),
          '_',
          substr(start = 7, stop = 7, Var2)
        )
      )  }

  if(what == "DBS"){

    a = a %>% dplyr::rename(Var1 = sbs, Var2 = variable) %>%
      mutate(
        substitution = paste0(substr(start = 1, stop = 2, Var2),">NN"),
        context = substr(start = 4, stop = 5, Var2)
      )
  }

  if(what == "ID"){

    a = a %>% dplyr::rename(Var1 = sbs, Var2 = variable) %>%
      mutate( Var2 = as.character(Var2),
              substitution = substr(start = 1, stop = nchar(Var2) - 2, Var2),
              context = substr(start = nchar(Var2), stop = nchar(Var2), Var2)
      )
  }

  if(what == "CNV"){

    a = a %>% dplyr::rename(Var1 = sbs, Var2 = variable) %>% rowwise() %>%
      mutate( Var2 = as.character(Var2),
              substitution =
                paste0(str_split(Var2,pattern = ":")[[1]][1],":",str_split(Var2,pattern = ":")[[1]][2]),
              context =  paste0(str_split(Var2,pattern = ":")[[1]][3])
      )
  }



  if(is.null(cls)){ cls = ggsci::pal_simpsons()(length(unique(a$Var1)))
  names(cls) = unique(a$Var1)
  }

  library(CNAqc)

  p = a  %>%
    ggplot() +
    geom_bar(aes(value, x = context, fill = Var1), stat = 'identity') +
    facet_grid(Var1 ~ factor(substitution,levels = a$substitution %>% unique()), scales = 'free') +
    CNAqc:::my_ggplot_theme() +
    scale_fill_manual(values = cls) +
    guides(fill = 'none')  +
    labs(y = "", title = "Signatures")

  if(!context) {
    p = p + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + labs(x = "")
  }

  p
}

plot_exposure_data = function(b,sample_name = T){

  cls = ggsci::pal_simpsons()(ncol(b))
  names(cls) = colnames(b)

  p =   ggplot(data = b %>% as.data.frame() %>% mutate(sample = rownames(b)) %>%
                 reshape2::melt() %>% dplyr::rename(Signature = variable),aes(x = sample, y  = value, fill = Signature)) +
    geom_bar(stat = "identity")  + ggplot2::scale_fill_manual(values = cls ) + labs(title = "Expsosure", y = "") +
    theme(axis.text.x = element_text(angle = 90))

  if (!sample_name) {
    p =  p +  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + labs(x = "")

  }

  p
}

get_similarity_scores = function(reconstr_data,real_data){

reconstr_data = reconstr_data[rownames(real_data),]
cosines = lapply(1:nrow(reconstr_data),function(i){

  tibble(id = rownames(real_data)[i], cos =  cosine.vector(reconstr_data[i,],real_data[i,]))

}) %>% bind_rows()

cosines
}


get_reconstructed_data = function(x){

 signatures = NULL
 if("catalogue_signatures" %in% names(x$fit)){ signatures = rbind(signatures,x$fit$catalogue_signatures) }
 if("denovo_signatures" %in% names(x$fit)){ signatures = rbind(signatures,x$fit$denovo_signatures)}

  signatures = signatures[colnames(x$fit$exposure),]

  as.matrix(x$fit$exposure[rownames(x$input$counts),]*rowSums(x$input$counts)) %*%  as.matrix(signatures) %>%  as.data.frame()

}


plot_group_data = function(data, what = "SBS", cls = c("Lung" = "forestgreen", 'Colorectal' = "purple3")){


a =  data %>% dplyr::mutate(sbs = rownames(data)) %>% as_tibble() %>%
    reshape2::melt()

  library(stringr)

  if(what == 'SBS'){
    a = a %>% dplyr::rename(Var1 = sbs, Var2 = variable) %>%
      mutate(
        substitution = paste0(substr(start = 3, stop = 3, Var2),">",substr(start = 5, stop = 5, Var2)),
        context = paste0(
          substr(start = 1, stop = 1, Var2),
          '_',
          substr(start = 7, stop = 7, Var2)
        )
      )  }

  if(what == "DBS"){

    a = a %>% dplyr::rename(Var1 = sbs, Var2 = variable) %>%
      mutate(
        substitution = paste0(substr(start = 1, stop = 2, Var2),">NN"),
        context = substr(start = 4, stop = 5, Var2)
      )
  }

  if(what == "ID"){

    a = a %>% dplyr::rename(Var1 = sbs, Var2 = variable) %>%
      mutate( Var2 = as.character(Var2),
              substitution = substr(start = 1, stop = nchar(Var2) - 2, Var2),
              context = substr(start = nchar(Var2), stop = nchar(Var2), Var2)
      )
  }

  if(what == "CNV"){

    a = a %>% dplyr::rename(Var1 = sbs, Var2 = variable) %>% rowwise() %>%
      mutate( Var2 = as.character(Var2),
              substitution =
                paste0(str_split(Var2,pattern = ":")[[1]][1],":",str_split(Var2,pattern = ":")[[1]][2]),
              context =  paste0(str_split(Var2,pattern = ":")[[1]][3])
      )
  }


  a %>%  dplyr::group_by(organ) %>%
    dplyr::mutate(n_muts=value/sum(value)) %>%
    ggplot() +
    geom_bar(aes(y = n_muts, x = context, fill = organ), stat = 'identity') +
    my_ggplot_theme()  +
    facet_grid(organ~ substitution, scales = 'free') + scale_fill_manual(values = cls) +
    theme(axis.ticks.x = element_blank(),axis.text.x = element_blank()
    ) + labs(x = '', y = "Mutation count")

}





# x = cnv fit

load("~/Documents/GitHub/simbasilica/processed_data/cna_counts.RData")
load("~/Documents/GitHub/simbasilica/processed_data/metadata.RData")
library(tidyr)

cnv_cosmic_catalogue = readRDS("~/Dropbox/Organoids_Accelerator/data/dbs-indel-cnv-Matrices/catalogues/CNA_cosmic_catalogue.rds")  %>%
  as.data.frame()

cna_data = cna_counts %>% t %>% as.data.frame() %>% mutate(PATIENT_ID = colnames(cna_counts)) %>%
  full_join(as.data.frame((metadata)), by = "PATIENT_ID") %>% filter(CANCER_TYPE %in% c("Breast Cancer","Glioblastoma" )) %>%
  drop_na()

rownames(cna_data) = cna_data$PATIENT_ID

b =   cna_data %>% filter(CANCER_TYPE == "Breast Cancer")
c =   cna_data %>% filter(CANCER_TYPE == "Glioblastoma" )

cna_data = rbind(b[1:50,],c[1:50,])

cancer_type = (cna_data$CANCER_TYPE == "Breast Cancer")*1

labels = cna_data %>% dplyr::select(c("PATIENT_ID","CANCER_TYPE")) %>% as_tibble() %>% rename(Sample = PATIENT_ID, groups = CANCER_TYPE)

cna_data = cna_data %>% dplyr::select(-c("PATIENT_ID","SAMPLE_ID","CANCER_TYPE")) %>% as.data.frame()



x = basilica::fit(x= cna_data, k=0:6, py=NULL,
                    reference_catalogue = cnv_cosmic_catalogue + 1e-18,
                    input_catalogue= cnv_cosmic_catalogue["CN17",] + 1e-18,lr = 0.01,steps = 500,groups = cancer_type,
                    reg_weight = 0)


saveRDS(x,"cnv_fit.rds")

 # plot sign
cls_data = ggsci::pal_nejm()(nrow(x$fit$catalogue_signatures))
names(cls_data) = rownames(x$fit$catalogue_signatures)

cnv_sign = plot_data_signatures(x$fit$catalogue_signatures,what = "CNV",context = F,cls = cls_data) + labs(title = 'Copy Number Signatures')
ggsave(cnv_sign,filename = "cnv_sign.png",height = 4,width = 7)

# plot exp

exp = plot_exposure(x,labels = labels) + labs(title = "Relative Exposure", x = "") +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())

ggsave(exp,filename = "cnv_exp.png",height = 4.5,width = 8.5)

# plot_exposure_data(x$fit$exposure,sample_name = F)
#
# # plot reconstr score
# cosines = get_similarity_scores(get_reconstructed_data(x),x$input$counts) %>%
#   full_join(metadata %>% as.data.frame() %>% dplyr::rename(id = PATIENT_ID,organ = CANCER_TYPE) %>% dplyr::select(id,organ), by = "id") %>%
#   drop_na()
#
#
# ggplot(cosines,aes(x = organ,y = cos, fill = organ)) + geom_boxplot() +
#   CNAqc:::my_ggplot_theme() +
#   labs(y = "Cosine Similarity", title = "Data reconstruction")  + theme(
#     axis.text.x = element_blank()
#   )
#
#
# ggplot(  x$fit$exposure %>% mutate(id = rownames(x$fit$exposure)) %>% full_join(metadata  %>%
#            as.data.frame() %>% dplyr::rename(id = PATIENT_ID,organ = CANCER_TYPE) , by = "id") %>% drop_na()  )
#
#
#
# # plot data
# rec_data = get_reconstructed_data(x)
#
# rec_data = rec_data %>% dplyr::mutate(id = rownames(rec_data))  %>%
#   full_join(metadata %>% as.data.frame() %>% dplyr::rename(id = PATIENT_ID,organ = CANCER_TYPE) %>% dplyr::select(id,organ),
#             by = "id") %>% drop_na()
#
# rownames(rec_data) = rec_data$id
#
# data = x$input$counts  %>% dplyr::mutate(id = rownames(x$input$counts))  %>%
#   full_join(metadata %>% as.data.frame() %>% dplyr::rename(id = PATIENT_ID,organ = CANCER_TYPE) %>% dplyr::select(id,organ),
#             by = "id") %>% drop_na()
#
# rownames(data) = data$id
#
# plot_group_data(rec_data, what = "CNV",cls = c("Non-Small Cell Lung Cancer" = "forestgreen", "Colorectal Cancer" = "purple3")) +
#   labs(title = "rec")
# plot_group_data(data, what = "CNV", cls = c("Non-Small Cell Lung Cancer" = "forestgreen", "Colorectal Cancer" = "purple3")) +
#   labs(title = "real")
#

# comparison plots real data

# lung and colorectal

x <- readRDS("~/Desktop/processed_data/test_fit.Rds")

exposure_lung = read.table("~/Desktop/processed_data/SBS_v2.03/organSpecificExposures/GEL/GEL-Lung_SBS_exposures_finalT.tsv",row.names = 1) %>%
  dplyr::select(-unassigned)
exposure_colorectal = read.table("~/Desktop/processed_data/SBS_v2.03/organSpecificExposures/GEL/GEL-Colorectal_SBS_exposures_finalT.tsv",
                                 row.names = 1)  %>% dplyr::select(-unassigned)

serena_signatures = read.table("~/Desktop/processed_data/SBS_v2.03/OrganSpecificSigs_GEL_SBS_v2.03.tsv",row.names = 1)

serena_signatures = serena_signatures[,c(colnames(exposure_lung),colnames(exposure_colorectal))] %>% t

exposure_serena = exposure_colorectal %>% mutate(id = rownames(exposure_colorectal),organ = "Colorectal") %>%
  full_join(exposure_lung %>% mutate(id = rownames(exposure_lung),organ = "Lung"), by = c("id","organ"))

exposure_serena[is.na(exposure_serena)] = 0

rownames(exposure_serena) =  exposure_serena$id

exposure_serena = exposure_serena %>% dplyr::select(-id)

exposure_serena = exposure_serena[rownames(x$input$counts),]


# data reconstruction

# organs = exposure_serena %>% dplyr::select(organ) %>% mutate(id = rownames(exposure_serena))
#
# exposure_serena =  exposure_serena %>% dplyr::select(-organ)
#
# serena_rec = as.matrix(exposure_serena*rowSums(x$input$counts)) %*%
#   as.matrix(serena_signatures[colnames(exposure_serena),]) %>%  as.data.frame() %>% mutate(id = rownames(exposure_serena)) %>%
#   full_join(organs, by = "id")
#
# basilica_rec = get_reconstructed_data(x) %>% mutate(id = rownames(x$input$counts)) %>% full_join(organs, by = "id")
#


# plot data

# data = x$input$counts %>% dplyr::mutate(id = rownames(x$input$counts)) %>% full_join(organs, by = "id")
# rownames(data) = data$id
#
# plot_group_data(serena_rec) + labs(title = "serena reconstruction")
# plot_group_data(basilica_rec) + labs(title = "basilica reconstruction")
# plot_group_data(data) + labs(title = "Real data")


# plot_group_exposure

exposure_serena %>% mutate(id = rownames(exposure_serena)) %>% full_join(organs, by = "id")





plot_exposure(x)



# plot cosines

# rownames(basilica_rec) = basilica_rec$id
# rownames(serena_rec) = serena_rec$id
#
# basilica_cosines = get_similarity_scores(basilica_rec %>% dplyr::select(-c("id","organ")),x$input$counts) %>% dplyr::rename(basilica = cos)
#
# serena_cosines = get_similarity_scores(serena_rec %>% dplyr::select(-c("id","organ")),x$input$counts) %>% dplyr::rename(serena = cos)
#
# cosines = full_join(basilica_cosines,serena_cosines, by = "id") %>% full_join(as_tibble(organs),by = "id")
#
# cosine_plot = ggplot(cosines %>% reshape2::melt(),aes(x = organ,y = value, fill = organ)) + geom_boxplot() +
#   CNAqc:::my_ggplot_theme() + facet_wrap(~variable) +
#   labs(y = "Cosine Similarity", title = "Data reconstruction")  + theme(
#     axis.text.x = element_blank()
#   )



# basilica signatures

basilica_sign = plot_signatures(x)

# serena signatures
# lung

plot_serena_sign = function(type,cls){

lung_sign = serena_signatures[grepl(rownames(serena_signatures),pattern = type),]

rownames(lung_sign) = stringr::str_remove(rownames(lung_sign),pattern = paste0("GEL.",type,"_common_"))
rownames(lung_sign) = stringr::str_remove(rownames(lung_sign),pattern = paste0("GEL.",type,"_rare_"))

names(cls) = rownames(lung_sign)
plot_data_signatures(lung_sign %>% as.data.frame(),cls = cls,context = F) + labs(title = paste0(type," Cancer"))

}

lung = plot_serena_sign("Lung", cls = c(ggsci::pal_simpsons()(16),ggsci::pal_jama()(4)))

col = plot_serena_sign("Colorectal",cls = c(ggsci::pal_simpsons()(16),ggsci::pal_nejm()(8),ggsci::pal_lancet()(4)) )






