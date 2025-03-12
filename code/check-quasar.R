
library(readr)
library(lme4)
library(dplyr)
library(tidyr)
library(MASS)

pheno_file <- "/home/jp2045/quasar-paper/data/nextflow/76/7ff46fa61e333fdda395a52f99e877/onek1k-Mono C-pheno-sum.txt"
cov_file <- "/home/jp2045/quasar-paper/data/nextflow/03/b84ae6e5fa73d6dce0cef6282407d3/onek1k-Mono C-all-covs.tsv"
grm_file <- "/home/jp2045/quasar-paper/data/nextflow/51/1f0fb4e777474356375d49d5c1ecb7/onek1k-grm.tsv"
pheno_data <- read_tsv(pheno_file)
cov_data <- read_tsv(cov_file)
grm_data <- read_tsv(grm_file)

gene_ind <- which(pheno_data$feature_id == "ENSG00000196199")

gene_data <- pheno_data[gene_ind, 2:ncol(pheno_data)] |>
  pivot_longer(
    cols = everything(),
    names_to = "sample_id",
    values_to = "count"
  )

all_data <- gene_data |>
  left_join(cov_data, by = "sample_id") |>
  arrange(sample_id)

fit_data <- all_data |>
  dplyr::select(-sample_id)

glm_nb_fit <- glm.nb(count ~ 0 + ., data = fit_data)
glm_fit <- glm(count ~ 0 + ., data = fit_data, family = "poisson")

# Fit GLMM.
sample_ind <- which(colnames(grm_data) %in% all_data$sample_id)
grm <- as.matrix(grm_data[sample_ind - 1, sample_ind])
colnames(grm) <- NULL

all_data

formula <- paste0(
  "count ~ sex + age + ",
  paste0("expr_pc", 1:5, collapse = " + "),
  " + ",
  paste0("geno_pc", 1:5, collapse = " + "),
  " + (1 | sample_id)"
) |> as.formula()

parsed_formula <- lFormula(
  formula = formula,
  data = all_data,
  control = lmerControl(check.nobs.vs.nlev = "ignore",
                        check.nobs.vs.nRE = "ignore"))

relmat <- list(sample_id = 2 * grm)
relfac <- relmat
flist <- parsed_formula$reTrms[["flist"]]
Ztlist <- parsed_formula$reTrms[["Ztlist"]]
stopifnot(all(names(relmat) %in% names(flist)))
asgn <- attr(flist, "assign")
for (i in seq_along(relmat)) {
  tn <- which(match(names(relmat)[i], names(flist)) == asgn)
  if (length(tn) > 1)
    stop("a relationship matrix must be associated",
         " with only one random effects term", call. = FALSE)
  relmat[[i]] <- Matrix(relmat[[i]], sparse = TRUE)
  relfac[[i]] <- chol(relmat[[i]])
  Ztlist[[i]] <- relfac[[i]] %*% Ztlist[[i]]
}
parsed_formula$reTrms[["Ztlist"]] <- Ztlist
parsed_formula$reTrms[["Zt"]] <- do.call(rbind, Ztlist)

parsed_formula$family <- poisson()
deviance_function <- do.call(mkGlmerDevfun, parsed_formula)
optimizer_output <- optimizeGlmer(deviance_function)

glmm_fit <- mkMerMod(
  rho = environment(deviance_function),
  opt = optimizer_output,
  reTrms = parsed_formula$reTrms,
  fr = parsed_formula$fr
)
glmm_fit
