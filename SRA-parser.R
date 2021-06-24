require(reutils)
require(XML)
require(dplyr)
require(stringr)

SRA_table <- function(GSE){

gse_death_parse <- function(GSE){
  db_list <- esearch(GSE,"gds") %>% 
    esummary(.) %>% content(., "parsed")
  gsm_db_df <- lapply(db_list, function(x) x$entryType == "GSM") %>%
    unlist() %>%
    db_list[.] %>%
    bind_rows() %>%
    select_if(function(x){!all(is.na(x))})
  names(gsm_db_df)[grepl("TargetObject",names(gsm_db_df))] <- "SRX"
  # combine all SRX
  Sys.sleep(2)
  l <- lapply(gsm_db_df$SRX, srx_death_parse)
  names(l) <- gsm_db_df$SRX
  df_temp <- bind_rows(l, .id = "SRX")
  full_join(gsm_db_df, df_temp)
}
srx_death_parse <- function(SRX){
  print(SRX)
  xml <- esearch(SRX,"sra") %>%
    efetch(.) %>% content()
  accs <- xpathApply(xml,"//RUN/@accession") %>% as.character()
  rl <- xpathApply(xml,"//Read[@count !=0]/@average") %>% as.numeric()
  len <- xpathApply(xml,"//Alternatives/@url") %>% as.character() %>% length()
  out <- data.frame(url = xpathApply(xml,"//Alternatives/@url") %>% as.character(),
                    org = xpathApply(xml,"//Alternatives/@org") %>% as.character()) %>%
    mutate(type = ifelse(grepl("fastq", url), "fastq", "sra")) %>%
    filter(type != "fastq")
  out$acc <- out$url %>% str_split(., "/") %>% lapply(., function(x) x[grepl("SRR", x)][1]) %>% unlist()
  out <- out %>% 
    group_by(org, acc) %>%
    slice(1) %>%
    ungroup() %>%
    dplyr::select(acc,org,url) %>%
    pivot_wider(names_from = org, values_from = url)
  lib <- xpathApply(xml,"//LIBRARY_LAYOUT")[[1]] %>%
    capture.output() %>% as.character() %>%
    .[2] %>% str_trim() %>%
    str_extract(., "[A-Z]+")
  out <- data.frame(acc = accs, read_length = rl) %>% left_join(out, .) %>%
    mutate(Library_type = lib) 
  Sys.sleep(2)
  out
}
out = gse_death_parse(GSE)
out = out %>% subset(., select = c(title, taxon, GCP, NCBI, read_length, Library_type )) %>%
  group_by(title,taxon, Library_type, read_length) %>%
  summarise(links = toString(GCP)) %>%
  ungroup() %>%
  mutate(links = gsub(links, pattern = ',', replacement = '')) %>%
  mutate(RSEM_ReadsPerGene = "NA",
         RSEM_SortedBAM = "NA",
         RSEM_BAMIndex = "NA",
         Kallisto_abundance = "NA",
         Kallisto_SortedBAM = "NA",
         Kallisto_BAMIndex = "NA")
colnames(out)[1] = 'entity:sample_id'

return(out)
}

# Example
out = SRA_table("GSE107218")
write_tsv(out, "out.tsv")
