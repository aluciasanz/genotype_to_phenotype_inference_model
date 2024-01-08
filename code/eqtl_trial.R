library(MatrixEQTL)

phage_nonsyn_table = phage_table[which(phage_table$nonsyn_mutation == 1), ]
rownames(phage_nonsyn_table) = 1:nrow(phage_nonsyn_table)

#### 4 types of phage table
phage_only_eqtl_m_allmut_allgen = phage_feature_binary_add_an_m_allmut_allgen
rownames(phage_only_eqtl_m_allmut_allgen) = c("phage_an", paste0("phage_mut_", 1:(nrow(phage_only_eqtl_m_allmut_allgen)-1)))

phage_only_eqtl_m_nonsyn_allgen = phage_feature_binary_add_an_m_nonsyn_allgen
rownames(phage_only_eqtl_m_nonsyn_allgen) = c("phage_an", paste0("phage_mut_nonsyn_", 1:(nrow(phage_only_eqtl_m_nonsyn_allgen)-1)))

phage_only_eqtl_m_allmut_uniqgen = phage_feature_binary_add_an_m_allmut_uniqgen
rownames(phage_only_eqtl_m_allmut_uniqgen) = c("phage_an", paste0("phage_mut_nonsyn_", 1:(nrow(phage_only_eqtl_m_allmut_uniqgen)-1)))

phage_only_eqtl_m_nonsyn_uniqgen = phage_feature_binary_add_an_m
rownames(phage_only_eqtl_m_nonsyn_uniqgen) = c("phage_an", paste0("phage_mut_nonsyn_", 1:(nrow(phage_only_eqtl_m_nonsyn_uniqgen)-1)))

#### 2 types of bac table
bac_only_eqtl_m_allmut_allgen = bac_feature_binary_add_an_m_allmut_allgen
rownames(bac_only_eqtl_m_allmut_allgen) = c("bac_an", paste0("bac_mut_", 1:(nrow(bac_only_eqtl_m_allmut_allgen)-1)))

bac_only_eqtl_m_allmut_uniqgen = bac_feature_binary_add_an_m
rownames(bac_only_eqtl_m_allmut_uniqgen) = c("bac_an", paste0("bac_mut_", 1:(nrow(bac_only_eqtl_m_allmut_uniqgen)-1)))


#### geno-table construction 1st order 2nd order

# 1. phage all mut all gen bac all mut all gen
p_aa_b_aa_1st <- construct_eqtl_feature_table_1st_inworking(phage_only_eqtl_m_allmut_allgen, bac_only_eqtl_m_allmut_allgen)
p_aa_b_aa_2nd <- construct_eqtl_feature_table_2nd(phage_only_eqtl_m_allmut_allgen, bac_only_eqtl_m_allmut_allgen)

# # 2. phage all mut uniq gen bac all mut all gen
# p_au_b_aa_1st <- p_na_b_aa_1st(phage_only_eqtl_m_allmut_uniqgen, bac_only_eqtl_m_allmut_allgen)
# p_au_b_aa_2nd <- construct_eqtl_feature_table_2nd(phage_only_eqtl_m_allmut_uniqgen, bac_only_eqtl_m_allmut_allgen)

# 3. phage nonsyn mut all gen bac all mut all gen
p_na_b_aa_1st <- construct_eqtl_feature_table_1st_inworking(phage_only_eqtl_m_nonsyn_allgen, bac_only_eqtl_m_allmut_allgen)
p_na_b_aa_2nd <- construct_eqtl_feature_table_2nd(phage_only_eqtl_m_nonsyn_allgen, bac_only_eqtl_m_allmut_allgen)
  
# # 4. phage nonsyn mut uniq gen bac all mut all gen
# p_nu_b_aa_1st <- p_na_b_aa_1st(phage_only_eqtl_m_nonsyn_uniqgen, bac_only_eqtl_m_allmut_allgen)
# p_nu_b_aa_2nd <- construct_eqtl_feature_table_2nd(phage_only_eqtl_m_nonsyn_uniqgen, bac_only_eqtl_m_allmut_allgen)

# # 5. phage all mut all gen bac all mut uniq gen
# p_aa_b_au_1st <- p_na_b_aa_1st(phage_only_eqtl_m_allmut_allgen, bac_only_eqtl_m_allmut_uniqgen)
# p_aa_b_au_2nd <- construct_eqtl_feature_table_2nd(phage_only_eqtl_m_allmut_allgen, bac_only_eqtl_m_allmut_uniqgen)

# 6. phage all mut uniq gen bac all mut uniq gen
p_au_b_au_1st <- construct_eqtl_feature_table_1st_inworking(phage_only_eqtl_m_allmut_uniqgen, bac_only_eqtl_m_allmut_uniqgen)
p_au_b_au_2nd <- construct_eqtl_feature_table_2nd(phage_only_eqtl_m_allmut_uniqgen, bac_only_eqtl_m_allmut_uniqgen)

# # 7. phage nonsyn mut all gen bac all mut uniq gen
# p_na_b_au_1st <- p_na_b_aa_1st(phage_only_eqtl_m_nonsyn_allgen, bac_only_eqtl_m_allmut_uniqgen)
# p_na_b_au_2nd <- construct_eqtl_feature_table_2nd(phage_only_eqtl_m_nonsyn_allgen, bac_only_eqtl_m_allmut_uniqgen)

# 8. phage nonsyn mut uniq gen bac all mut uniq gen
p_nu_b_au_1st <- construct_eqtl_feature_table_1st_inworking(phage_only_eqtl_m_nonsyn_uniqgen, bac_only_eqtl_m_allmut_uniqgen)
p_nu_b_au_2nd <- construct_eqtl_feature_table_2nd(phage_only_eqtl_m_nonsyn_uniqgen, bac_only_eqtl_m_allmut_uniqgen)



#### pheno-table and label construction 

p_allgen_name_order = colnames(phage_only_eqtl_m_allmut_allgen)
p_uniqgen_name_order = colnames(phage_only_eqtl_m_allmut_uniqgen)
b_allgen_name_order = colnames(bac_only_eqtl_m_allmut_allgen)
b_uniqgen_name_order = colnames(bac_only_eqtl_m_allmut_uniqgen)

p_allgen_b_allgen_pheno_m = pheno_matrix_add_an[, p_allgen_name_order]
p_allgen_b_allgen_pheno_m = p_allgen_b_allgen_pheno_m[b_allgen_name_order, ]
p_uniqgen_b_uniqgen_pheno_m = pheno_matrix_add_an_ordered_uniq_avg[, p_uniqgen_name_order]
p_uniqgen_b_uniqgen_pheno_m = p_uniqgen_b_uniqgen_pheno_m[b_uniqgen_name_order, ]

p_allgen_b_allgen_pheno_name = as.vector(as.matrix(sapply(colnames(p_allgen_b_allgen_pheno_m), function(x) paste0(x, ":", rownames(p_allgen_b_allgen_pheno_m)))))
p_uniqgen_b_uniqgen_pheno_name = as.vector(as.matrix(sapply(colnames(p_uniqgen_b_uniqgen_pheno_m), function(x) paste0(x, ":", rownames(p_uniqgen_b_uniqgen_pheno_m)))))

p_allgen_b_allgen_pheno_l = as.vector(as.matrix(p_allgen_b_allgen_pheno_m))
names(p_allgen_b_allgen_pheno_l) = p_allgen_b_allgen_pheno_name
p_uniqgen_b_uniqgen_pheno_l = as.vector(as.matrix(p_uniqgen_b_uniqgen_pheno_m))
names(p_uniqgen_b_uniqgen_pheno_l) = p_uniqgen_b_uniqgen_pheno_name

label_allgen = p_allgen_b_allgen_pheno_l
label_uniqgen = p_uniqgen_b_uniqgen_pheno_l


p_phen_allgen = as.matrix(p_allgen_b_allgen_pheno_m)
b_phen_allgen = t(as.matrix(p_allgen_b_allgen_pheno_m))

p_phen_uniqgen = as.matrix(p_uniqgen_b_uniqgen_pheno_m)
b_phen_uniqgen = t(as.matrix(p_uniqgen_b_uniqgen_pheno_m))
