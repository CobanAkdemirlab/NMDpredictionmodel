# ==============================================================================
# NMDpredictionmodel - upstream R dependencies
#
# Used by:
#   Dataset_extraction/
#   Annotation/
#   Features/
#
# Run with:
#
#   Rscript dependencies/install_R_dependencies.R
#
# ==============================================================================


# ------------------------------------------------------------------------------
# CRAN packages
# ------------------------------------------------------------------------------

cran_packages <- c(
    "dplyr",
    "tidyr",
    "data.table",
    "readr",
    "stringr",
    "tibble",
    "purrr",
    "remotes"
)


installed <- rownames(
    installed.packages()
)


for (pkg in cran_packages) {

    if (!pkg %in% installed) {

        message(
            "Installing CRAN package: ",
            pkg
        )

        install.packages(
            pkg,
            repos = "https://cloud.r-project.org"
        )
    }
}


# ------------------------------------------------------------------------------
# Bioconductor manager
# ------------------------------------------------------------------------------

if (!requireNamespace(
    "BiocManager",
    quietly = TRUE
)) {

    install.packages(
        "BiocManager",
        repos = "https://cloud.r-project.org"
    )
}


# ------------------------------------------------------------------------------
# Bioconductor packages
#
# Used for VCF handling, genomic ranges, transcript structures,
# genome sequence manipulation, and GENCODE annotation.
# ------------------------------------------------------------------------------

bioc_packages <- c(
    "S4Vectors",
    "IRanges",
    "GenomicRanges",
    "GenomeInfoDb",
    "GenomicFeatures",
    "Biostrings",
    "VariantAnnotation",
    "rtracklayer"
)


for (pkg in bioc_packages) {

    if (!requireNamespace(
        pkg,
        quietly = TRUE
    )) {

        message(
            "Installing Bioconductor package: ",
            pkg
        )

        BiocManager::install(
            pkg,
            ask = FALSE,
            update = FALSE
        )
    }
}


# ------------------------------------------------------------------------------
# aenmd
#
# Used during gnomAD, ClinVar, and GREGoR extraction to identify
# variants predicted to create premature termination codons.
# ------------------------------------------------------------------------------


if (!requireNamespace(
    "aenmd.data.ensdb.v105",
    quietly = TRUE
)) {

    message(
        "Installing aenmd annotation data..."
    )

    remotes::install_github(
        repo = "kostkalab/aenmd_data",
        subdir = "aenmd.data.ensdb.v105",
        upgrade = "never"
    )
}


# ------------------------------------------------------------------------------
# aenmd
# ------------------------------------------------------------------------------

if (!requireNamespace(
    "aenmd",
    quietly = TRUE
)) {

    message(
        "Installing aenmd..."
    )

    remotes::install_github(
        "kostkalab/aenmd",
        dependencies = TRUE,
        upgrade = "never"
    )
}

# ------------------------------------------------------------------------------
# Verify installation
# ------------------------------------------------------------------------------

required_packages <- c(
    "aenmd.data.ensdb.v105",
    "aenmd",
    "S4Vectors",
    "IRanges",
    "GenomicRanges",
    "GenomeInfoDb",
    "GenomicFeatures",
    "Biostrings",
    "VariantAnnotation",
    "rtracklayer",
    "dplyr",
    "tidyr",
    "data.table"
)


message("")
message("============================================================")
message("R DEPENDENCY CHECK")
message("============================================================")


failed <- character()


for (pkg in required_packages) {

    if (requireNamespace(
        pkg,
        quietly = TRUE
    )) {

        message(
            "OK     ",
            pkg,
            " ",
            as.character(
                packageVersion(pkg)
            )
        )

    } else {

        message(
            "FAILED ",
            pkg
        )

        failed <- c(
            failed,
            pkg
        )
    }
}


if (length(failed) > 0) {

    stop(
        "Missing R packages: ",
        paste(
            failed,
            collapse = ", "
        )
    )
}


message("")
message(
    "All upstream R dependencies are available."
)
