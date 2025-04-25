#' @title WHPT-ASPT Score Calculator
#' @description Returns the WHPT-ASPT score for each sample submitted.
#' @importFrom dplyr filter mutate group_by summarize select left_join case_when
#' @importFrom tidyr pivot_longer
#' @param x Benthic Macroinvertebrate dataframe (csv or data.frame).
#' @param Sample_ID Name of column containing Sample Name.
#' @details
#' The function uses reference tables to score the abundance and presence-absence
#' data, and combines these into an overall score for each sample. A template
#' for the data input can be found in the Google Sheets template here:
#' \url{https://docs.google.com/spreadsheets/d/1Z0Y5O0bQXqJ_14iYOwn4UVe4L4MdhdfV7HSwolNhvdg/edit?usp=sharing}.
#' @examples
#' # Example Usage:
#' # Load the example dataset included in the package
#' data("test_whpt_dataset")
#' WHPT_ASPT_Score(test_whpt_dataset)
#' @export
WHPT_ASPT_Score <- function(x, Sample_Label = NULL) {
  data("Ref_Long", envir = environment())
  data("WHPT_Families", envir = environment())
  data("P_A_WHPT_Scores", envir = environment())
  data("test_whpt_dataset", envir = environment())

  # Determine sample ID column
  if (!is.null(Sample_Label)) {
    if (!(Sample_Label %in% colnames(x))) {
      stop("The specified 'Sample_Label' does not exist in the dataset.")
    }
    message("Using column: ", Sample_Label)
  } else {
    message("No 'Sample_Label' specified. Using the first column by default.")
    Sample_Label <- colnames(x)[1]
  }

  # Rename Sample_Label column to Site_Label
  colnames(x)[colnames(x) == Sample_Label] <- "Site_Label"

  # Pivot to long format
  data_long <- tidyr::pivot_longer(
    data = x,
    cols = -Site_Label,
    names_to = "Family",
    values_to = "Abundance"
  )

  # Filter out zero abundances
  data_long <- dplyr::filter(data_long, Abundance != 0)

  # Abundance-based scoring --------------------------------------------
  abundance_data <- data_long %>%
    dplyr::mutate(
      A_Cat = dplyr::case_when(
        Abundance >= 1 & Abundance <= 9 ~ "A1",
        Abundance >= 10 & Abundance <= 99 ~ "A2",
        Abundance >= 100 & Abundance <= 999 ~ "A3",
        Abundance >= 1000 & Abundance <= 10000 ~ "A4",
        TRUE ~ "Other"
      )
    )

  # Check for unmatched families
  invalid_families <- setdiff(unique(abundance_data$Family), unique(Ref_Long$Family))
  if (length(invalid_families) > 0) {
    warning("The following families are not found in the reference database:\n",
            paste(invalid_families, collapse = ", "))
    return(invalid_families)
  } else {
    message("All Abundance-weighted families found in the reference database!")
  }

  # Join with reference table
  abundance_joined <- dplyr::left_join(
    dplyr::select(abundance_data, -Abundance),
    Ref_Long,
    by = c("Family", "A_Cat")
  )

  if (any(is.na(abundance_joined$Abundance_Value))) {
    warning("Some Family/A_Cat combinations did not match in the reference table.")
  }

  abundance_score <- abundance_joined %>%
    dplyr::group_by(Site_Label) %>%
    dplyr::summarize(ASPT = mean(Abundance_Value, na.rm = TRUE))

  # Presence-absence scoring --------------------------------------------
  pa_data <- data_long %>%
    dplyr::distinct(Site_Label, Family)

  # Check for unmatched families in presence-absence data
  invalid_families_pa <- setdiff(unique(pa_data$Family), unique(P_A_WHPT_Scores$Family))
  if (length(invalid_families_pa) > 0) {
    warning("The following families in presence-absence data are not found in the PA reference database:\n",
            paste(invalid_families_pa, collapse = ", "))
    return(invalid_families_pa)
  } else {
    message("All PA families found in the presence-absence reference database!")
  }



  PA_joined <- dplyr::left_join(
    pa_data,
    P_A_WHPT_Scores,
    by = c("Family")
  )

  PA_score <- PA_joined %>%
    dplyr::group_by(Site_Label) %>%
    dplyr::summarize(PA = mean(P_A_Value, na.rm = TRUE))


  # Combine both results --------------------------------------------
  final_scores <- dplyr::left_join(abundance_score, PA_score, by = "Site_Label")

  return(final_scores)
}
