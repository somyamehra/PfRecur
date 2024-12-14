library(readxl)
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
library(reshape2)
library(tidyr)
library(stringr)
library(readr)

setwd(this.path::here())

# ============================ HELPER FUNCTIONS ============================

add_prefix <- function(allele, prefix) {
  ifelse(!is.na(allele), gsub(", ", paste0("/", prefix), paste0(prefix, allele)), NA)
}

calc_MOI <- function(allele) {
  ifelse(allele!="", str_count(allele, "/")+1, 0)
}

# ============================ PARSE GENOTYPES ============================

Angola_TES <- bind_rows(read_xlsx("Dimbu_Angola_TES_genotypes.xlsx", skip=3),
                         read_xlsx("Dimbu_Angola_TES_genotypes.xlsx", sheet=2, skip=3)) %>%
  as.data.frame()

MARKERS <- c(M313="M313", M383="M383", TA1="TA1", POLYA="POLYA",
             PFPK2="PFPK2", M2490="M2490", TA109="TA109") 

STUDY_SITES <- c("Benguela", "Lunda Sul", "Zaire")

Angola_TES_data <- Angola_TES %>% rename(Sample_ID=`Sample.ID`) %>%
  extract(Sample_ID, c("Code", "Timepoint"), "([A-Z]{2}[0-9-]{6,})(D[0-9]+)", remove = F) %>%
  unite(TA109, colnames(Angola_TES)[grepl("TA109", colnames(Angola_TES))], sep=", ", na.rm=T) %>%
  unite(M313, colnames(Angola_TES)[grepl("313", colnames(Angola_TES))], sep=", ", na.rm=T) %>%
  unite(M383, colnames(Angola_TES)[grepl("383", colnames(Angola_TES))], sep=", ", na.rm=T) %>%
  unite(TA1, colnames(Angola_TES)[grepl("TA1_", colnames(Angola_TES))], sep=", ", na.rm=T) %>%
  unite(POLYA, colnames(Angola_TES)[grepl("POLYA", colnames(Angola_TES))], sep=", ", na.rm=T) %>%
  unite(PFPK2, colnames(Angola_TES)[grepl("PFPK2", colnames(Angola_TES))], sep=", ", na.rm=T) %>%
  unite(M2490, colnames(Angola_TES)[grepl("2490", colnames(Angola_TES))], sep=", ", na.rm=T) %>%
  mutate(TA109=add_prefix(TA109, ""),
         M313=add_prefix(M313, ""),
         M383=add_prefix(M383, ""),
         TA1=add_prefix(TA1, ""),
         POLYA=add_prefix(POLYA, ""),
         PFPK2=add_prefix(PFPK2, ""),
         M2490=add_prefix(M2490, "")) %>%
  mutate(MOI=pmax(calc_MOI(TA109), calc_MOI(M313), calc_MOI(M383), calc_MOI(TA1),
                  calc_MOI(POLYA), calc_MOI(PFPK2), calc_MOI(M2490))) %>%
  subset(MOI>0)
rownames(Angola_TES_data) <- Angola_TES_data$Sample_ID

Angola_TES_pairs <- Angola_TES_data %>% 
  mutate(Timepoint=ifelse(grepl("D0", Timepoint), "D_0", "D_REC")) %>%
  reshape2::dcast(formula=Code~Timepoint, value.var="Sample_ID") %>%
  mutate(Dimbu_posterior=Angola_TES_data[D_0, "Plucinski_posterior"])

Angola_TES_pairs <- Angola_TES_pairs %>% subset(!is.na(D_REC))

Angola_TES_marker_set <- lapply(MARKERS, function(m) {
  sort(unique(setdiff(unlist(strsplit(Angola_TES_data[,m], "/")), "NA")))})

Angola_TES_genotype_matrix <- lapply(MARKERS, function(m) {
  keep_samples <- subset(Angola_TES_data[, "Sample_ID"], 
                         Angola_TES_data[, m]!="")
  GT <- matrix(0, nrow=length(keep_samples), ncol=length(Angola_TES_marker_set[[m]]),
               dimnames = list(keep_samples, Angola_TES_marker_set[[m]]))
  for (indiv in keep_samples) {
    GT[indiv, unlist(setdiff(strsplit(Angola_TES_data[indiv, m], "/"), "NA"))] <- 1
  }
  return(GT)
})

# ============================ PARSE ISOLATES ============================

day_0_isolates <- 
  split(subset(Angola_TES_data$Sample_ID, Angola_TES_data$Timepoint=="D0"), 
        subset(substr(Angola_TES_data$Site, 1, 1), Angola_TES_data$Timepoint=="D0"))

Angola_TES_isolates <- list()
for (i in 1:nrow(Angola_TES_pairs)) {
  Angola_TES_isolates[[i]] <- list(recurrent=Angola_TES_pairs[i, "D_REC"],
                                   ref_C=Angola_TES_pairs[i, "D_0"],
                                   ref_I=setdiff(day_0_isolates[[substr(Angola_TES_pairs[i, "D_0"], 1, 1)]], 
                                                 Angola_TES_pairs[i, "D_0"]))
}

usethis::use_data(Angola_TES_genotype_matrix, Angola_TES_isolates,
                  internal = FALSE, overwrite = TRUE)



