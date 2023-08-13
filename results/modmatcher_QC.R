library(Modmatcher)

cwd <- getwd()
cohorts <- list('COAD')#, 'KIRC', 'KIRP', 'LUSC', 'PCPG', 'READ', 'STAD')

for (cohort in cohorts) {
    expr_data <- read.csv(paste0(cwd,"/datasets/", paste0("TCGA-", cohort, ".htseq_fpkm.tsv")), header = TRUE, row.names=1, sep = "\t")
    pheno_data <- read.csv(paste0(cwd,"/datasets/", paste0("TCGA-", cohort, ".GDC_phenotype.tsv")), header = TRUE, row.names=1, sep = "\t")#, fill=TRUE)

    # Remove everything after the first dot (".") in gene names and the dot itself
    first_row <- colnames(expr_data)
    cleaned_gene_names <- sub("\\..*", "", first_row)
    colnames(expr_data) <- cleaned_gene_names

    # Select specific gene columns
    selected_genes <- c("ENSG00000129824", "ENSG00000067048")
    expr_data_selected <- expr_data[selected_genes]
    #expr_data_selected <- expr_data_selected[expr_data_selected["ENSG00000129824"] != 0,]
    #expr_data_selected <- expr_data_selected[expr_data_selected["ENSG00000067048"] != 0,]

    # Select phenotype columns containing 'Sex', 'sex', or 'gender'
    pheno_columns <- colnames(pheno_data)
    selected_pheno_columns <- pheno_columns[grep("Sex|gender", pheno_columns, ignore.case=TRUE)]
    pheno_data_selected <- pheno_data[c(selected_pheno_columns)]

    # Merge the two tables by sample identifier
    merged_data <- merge(expr_data_selected, pheno_data_selected, by='row.names')
    merged_data <- merged_data[!is.na(merged_data[,4]),]
    row.names(merged_data) <- merged_data[,1]
    merged_data <- merged_data[,-1]

    # Extract gender information from the merged data
    ann_gender <- merged_data[,3]
    gender <- ifelse(is.na(ann_gender),3,ifelse(ann_gender=="Female" | ann_gender=='female',1,2))

    # Perform analysis using gender information
    gen_tab <- gender_sig(t(merged_data[,1:2]), gender)
    print(gen_tab)
    write.csv(gen_tab, paste0(cohort, '_modmatcher_gen_tab.csv'), row.names=FALSE)

    # Extract the most separating and second best separating genes
    first <- 'ENSG00000129824'  # The most separating genes between genders: RPS4Y1
    second <- 'ENSG00000067048'  # The second best separating genes between genders: DDX3Y

    # Check gender classification/separation
    gender_check(t(merged_data[,1:2]), first, second, gender)
} 