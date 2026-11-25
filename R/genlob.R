# Thinning response function
thin_response <- function(age, thin_age, thin_ratio, r, k) {
  response <- ifelse(age - thin_age > k | age - thin_age < 0, 1,
                     (1/thin_ratio)^((r*(k*(age-thin_age) - (age-thin_age)^2)) / age^2))

  # prevents accidental triggering
  if(thin_age < 0 || is.na(thin_age)) {response <- 1}

  return(response)
}

# Fertilization response function
fert_response <- function(age, fert_age, N, P, delta0, delta1, beta, gamma = 1) {
  yst <- age - fert_age
  response <- ifelse(yst <= 0, 1,
                     1 + ((delta0 + delta1 * P) * N) * (gamma / beta) *
                       (yst / beta)^(gamma - 1) * exp(-(yst / beta)^gamma)
  )

  # prevents accidental triggering
  if(fert_age < 0 || is.na(fert_age)) {response <- 1}

  return(response)
}

#' Calculates trees per hectare, basal area, and height from age1 to age2
#'
#' Takes into account growth response to fertilization and thinning.
#' Can run multiple scenarios by providing a vector of inputs,
#' but all vectors must be the same length or length 1.
#'
#' Model from paper by Gyawali and Burkhart, 2015.
#' Trained on loblolly pine plantation data in the southeastern US, ages 8 to 25.
#'
#'
#' @param age1 A numeric vector of positive integers with 1+ elements.
#'    Starting prediction age.
#' @param age2 A numeric vector of positive integers with 1+ elements.
#'    Ending prediction age.
#' @param n A numeric vector of positive numbers with 1+ elements.
#'    Number of trees per hectare.
#' @param ht A numeric vector of positive numbers with 1+ elements.
#'    Average dominant height in meters.
#' @param ba A numeric vector of positive numbers with 1+ elements.
#'    Average basal area in meters squared per hectare.
#' @param site_index A numeric vector of positive numbers with 1+ elements.
#'    Average dominant height in meters at a base age of 25 years.
#' @param thin_age A numeric vector of positive integers with 1+ elements.
#'    Year that the stand is thinned at.
#' @param thin_ratio A numeric vector of numbers between 0 and 1 with 1+ elements.
#'    Basal area before thinning divided by basal area after thinning.
#' @param fert_age A numeric vector of positive integers with 1+ elements.
#'    Year that the stand is fertilized at.
#' @param kgNperha A numeric vector of positive numbers with 1+ elements.
#'    Kilograms of nitrogen applied per hectare.
#' @param P A logical vector with 1 or more elements.
#'    If phosphorous was applied or not.
#' @param units A character vector with 1+ elements. 'metric' or 'imperial'
#' @param baseline_type A character vector with 1+ elements.
#'    'thin', 'fert' or 'base' chooses which study's baseline parameters to use.
#' @param BA_base_params A numeric vector with 2 elements.
#'    Default: 0.1488, 0.2978
#' @param BA_thin_params A numeric vector with 2 elements.
#'    Default: 0.1190, 30.8038
#' @param BA_fert_params A numeric vector with 3 elements.
#'    Default: 0.000738, 0.000538, 6.9419
#' @param H_base_params A numeric vector with 2 elements.
#'    Default: 0.0521, 0.8313
#' @param H_thin_params A numeric vector with 2 elements.
#'    Default: -0.1615, 5.4613
#' @param H_fert_params A numeric vector with 4 elements.
#'    Default: 0.000256, 0.000192, 7.7046, 2.9453
#' @param survival_params A numeric vector with 4 elements.
#'    Default: -1.1955, 1.22e-8, 1.8965, 3.2693
#'
#' @returns A dataframe with columns runID, age, n, ht, and ba.
#' Multiple scenarios appended to bottom of dataframe with different runIDs.
#' @export
#'
#' @examples
#' genlob(age1 = 10, age2 = 30, n = 1000, ht = 10, ba = 10, site_index = 18)
#'
#' genlob(age1 = 10, age2 = 30, n = 1000, ht = 10, ba = 10, site_index = 18,
#' thin_age = c(NA, NA, 15, 15), thin_ratio = c(1, 1, 0.5, 0.4),
#' fert_age = c(15, 15, NA, 15), kgNperha = c(112, 224, 0, 336), P = F)
#'
genlob <- function(
    # User inputs
  age1, age2, n, ht = NA, ba, # ages, number of trees, height, and ba required
  site_index = NA, # can use site_index or ht
  thin_age = NA, thin_ratio = 1, # if a thin takes place
  fert_age = NA, kgNperha = 0, P = FALSE, # if a fertilization takes place.
  units = "metric", # metric or imperial conversion to add in later
  baseline_type = "thin", # type of baseline parameter. "thin", "fert", "base"

  # Model parameters
  BA_base_params = NA, BA_thin_params = NA, BA_fert_params = NA,
  H_base_params = NA, H_thin_params = NA, H_fert_params = NA,
  survival_params = NA) {

  # Test non optional inputs so list doesn't break
  if(missing(age1) | missing(age2) | missing(n) | missing(ba) ||
     anyNA(age1) || anyNA(age2) || anyNA(n) || anyNA(ba)) {
    stop("The inputs 'age1', 'age2', 'n', and 'ba' MUST be supplied and not NA")
  }

  inputs <- list(age1 = age1, age2 = age2, n = n, ht = ht, ba = ba,
                 site_index = site_index,
                 thin_age = thin_age, thin_ratio = thin_ratio,
                 fert_age = fert_age, kgNperha = kgNperha, P = P,
                 units = units)

  lengths <- lengths(inputs)
  num_runs <- max(lengths)
  wrong <- names(inputs[lengths != 1 & lengths!= num_runs])
  if(length(wrong) > 0) {
    stop("All inputs must be either length 1 or the same max length of ",
         max(lengths)," \n The following inputs do not match this pattern: ",
         paste("\n  - ", wrong), call. = F)
  }

  # Convert to dataframe so all rows equal
  inputs <- data.frame(inputs)
  inputs$runID <- (1:num_runs)

  with(inputs, {
    if(any(age1 > age2))
      stop("'age1' must be less than or equal to 'age2'", call. = F)
    if(any(age1 <= 0))
      stop("'age1' must be greater than 0", call. = F)
    if(any(n <= 0))
      stop("'n' must be greater than 0", call. = F)
    if(any(!is.na(ht) & ht <= 0))
      stop("'ht' must be greater than 0", call. = F)
    if(any(ba <= 0))
      stop("'ba' must be greater than 0", call. = F)
    if(any(is.na(ht) & is.na(site_index)))
      stop("Must provide at least 1 of 'site_index' or 'ht'", call. = F)
    if(any(!is.na(site_index) & site_index <= 0))
      stop("'site_index' must be greater than 0", call. = F)
    if(any(is.na(thin_ratio) | thin_ratio > 1 | thin_ratio <= 0))
      stop("'thin_ratio' must fall between 0 and 1 or be equal to 1", call. = F)
    if(any(!is.na(thin_age) & (thin_age <= age1 | thin_age >= age2)))
      stop("'thin_age' must fall in between or be equal to 'age1' and 'age2'", call. = F)
    if(any(!is.na(thin_age) & thin_ratio == 1))
      warning("A 'thin_age' is provided but 'thin_ratio' is still 1", call. = F)
    if(any(thin_ratio != 1 & is.na(thin_age)))
      warning("A 'thin_ratio' is provided but thin_age is NA. No effect will occur", call. = F)
    if(any(is.na(kgNperha) | kgNperha < 0))
      stop("'kgNperha' must be greater than or equal to 0", call. = F)
    if(any(is.na(P) | !is.logical(P)))
      stop("'P' must be TRUE or FALSE", call. = F)
    if(any(!is.na(fert_age) & (fert_age <= age1 | fert_age >= age2)))
      stop("'fert_age' must fall in between or be equal to 'age1' and 'age2'", call. = F)
    if(any(!is.na(fert_age) & kgNperha == 0))
      warning("A 'fert_age' is provided but 'kgNperha' is still 0", call. = F)
    if(any((kgNperha != 0 | P != F) & is.na(fert_age)))
      warning("'kgNperha' and/or 'P' is provided but fert_age is NA. No effect will occur", call. = F)
    if(any(!units %in% c("metric", "imperial")))
      stop("'units' must be one of 'metric' or 'imperial'", call. = F)
  })

  # Ensure baseline_type is length 1 to prevent errors
  if(length(baseline_type) > 1)
    stop("'baseline_type' must be length 1. \n Calling different model ",
         "parameters in one run is not currently supported", call. = F)

  # Make sure they define the correct baseline type if modifying baseline parameters
  if((!missing(BA_base_params) | !missing(H_base_params) | !missing(survival_params)) &
     tolower(baseline_type) != "custom") {
    warning(paste0("You have defined custom baseline parameters but not selected 'baseline_type' custom. \n",
                   "Overriding with the parameters from type = ", baseline_type))
  }

  # Helper function to check that parameter arguments aren't missing or invalid
  is.valid <- function(vector, length) {
    if(missing(vector)) return(FALSE)
    is.numeric(vector) & length(vector) == length
  }

  # Define baseline parameters
  if (tolower(baseline_type) == "thin") {
    BA_base_params <- c(0.1488, 0.2978)
    H_base_params <- c(0.0521, 0.8313)
    survival_params <- c(-1.1955, 1.22e-8, 1.8965, 3.2693)

  } else if (tolower(baseline_type) == "fert") {
    BA_base_params <- c(0.1263, 0.2998)
    H_base_params <- c(0.0987, 0.4969)
    # defaults back to 'thin' for survival parameters
    survival_params <- c(-1.1955, 1.22e-8, 1.8965, 3.2693)
    message("Survival parameters were not estimated from 'fert' data. \n",
            "Defaulting back to baseline 'thin' for survival parameters only.")

  } else if (tolower(baseline_type) == "base") {
    BA_base_params <- c(0.1344, 0.3324)
    H_base_params <- c(0.0365, 0.8997)
    survival_params <- c(-1.0023, 5.573e-8, 1.5661, 3.3207)

  } else if (tolower(baseline_type) == "custom") {
    # Check is supplied custom values are valid
    if(!is.valid(BA_base_params, 2) | !is.valid(H_base_params, 2) |
       !is.valid(survival_params, 4)) {
      stop("baseline_type 'custom' selected with invalid arguments. \n",
           "You must provide a valid argument for: \n",
           "'BA_base_params', 'H_base_params', and 'survival_params'. \n",
           "These should be vectors of 2, 2, and 4 numbers, respectively")
    } else { # If valid, let users know what they input
      message(paste0("baseline_type 'custom' selected. Using the following parameters:",
                     "\nBA_base_params = ", toString(BA_base_params),
                     "\nH_base_params = ", toString(H_base_params),
                     "\nsurvival_params = ", toString(survival_params)))
    }
  } else {
    stop("baseline_type is not recognized. Must be one of: 'thin', 'fert', 'base', 'custom'")
  }

  if(missing(BA_thin_params)) {
    BA_thin_params <- c(0.1190, 30.8038)
  } else if (is.valid(BA_thin_params, 2)) {
    message("BA_thin_params modified to ", toString(BA_thin_params))
  } else {
    stop("Invalid 'BA_thin_params' argument. \n",
         "'BA_thin_params' needs to be a vector of 2 numbers", call. = F)
  }

  if(missing(BA_fert_params)) {
    BA_fert_params <- c(0.000738, 0.000538, 6.9419)
  } else if (is.valid(BA_fert_params, 3)) {
    message("BA_fert_params modified to ", toString(BA_fert_params))
  } else {
    stop("Invalid 'BA_fert_params' argument. \n",
         "'BA_fert_params' needs to be a vector of 3 numbers", call. = F)
  }

  if(missing(H_thin_params)) {
    H_thin_params <- c(-0.1615, 5.4613)
  } else if (is.valid(H_thin_params, 2)) {
    message("H_thin_params modified to ", toString(H_thin_params))
  } else {
    stop("Invalid 'H_thin_params' argument. \n",
         "'H_thin_params' needs to be a vector of 2 numbers", call. = F)
  }

  if(missing(H_fert_params)) {
    H_fert_params <- c(0.000256, 0.000192, 7.7046, 2.9453)
  } else if (is.valid(H_fert_params, 4)) {
    message("H_fert_params modified to ", toString(H_fert_params))
  } else {
    stop("Invalid 'H_fert_params' argument. \n",
         "'H_fert_params' needs to be a vector of 4 numbers", call. = F)
  }

  predict_ht <- function(age1, age2, ht) {
    ht * ((1 - exp(-H_base_params[1] * age2)) /
            (1 - exp(-H_base_params[1] * age1)))^(1 / H_base_params[2])
  }
  predict_ba <- function(age1, age2, ba) {
    ba * ((1 - exp(-BA_base_params[1] * age2)) /
            (1 - exp(-BA_base_params[1] * age1)))^(1 / BA_base_params[2])
  }

  # Store unit conversion variables
  acre_to_hectare <- 0.40468564224
  ft_to_m <- 0.3048
  m2_ha_to_ft2_ac <- 4.356

  combined_data <- list()

  for (r in inputs$runID) {
    single_run_data <- inputs[inputs$runID == r, ]
    combined_data[[r]] <- with(single_run_data, {

      initial_si <- site_index # save for reference in case it gets adjusted

      # save height_units for message output
      if(units == "imperial") {height_units <- "feet"}
      else {height_units <- "meters"}

      # predict site index if only height given
      if(is.na(site_index)) {
        site_index <- predict_ht(age1, 25, ht) # Doesn't include thin / fert effects
        message(paste0("Site index not provided for run ", r,
                       "\nEstimated for base age 25 using height of ",
                       ht, " ", height_units," and age ", age1,
                       ".\nResulting site index: (", height_units, "): ",
                       round(site_index, 4), "\n"))
      }

      # predict height if only site index given
      if(is.na(ht)) {
        ht <- predict_ht(25, age1, initial_si)
        message(paste0("Height not provided for run ", r,
                       "\nEstimated using provided site index of ",
                       initial_si, " ", height_units, " at base age 25",
                       ".\nResulting height (", height_units, "): ",
                       round(ht, 4), "\n"))
      }

      # Because height prediction was not unit dependent this can happen afterwards
      if(units == "imperial") {
        n <- n / acre_to_hectare
        ht <- ht * ft_to_m
        ba <- ba / m2_ha_to_ft2_ac
        site_index <- site_index * ft_to_m
      }

      valid_thin <- !is.na(thin_age) & thin_ratio != 1 & thin_age <= age2 & thin_age >= age1
      nrows <- length(age1:age2) + valid_thin
      # Set up vectors before loop

      age_vector <- numeric(nrows)
      n_vector <- numeric(nrows)
      ht_vector <- numeric(nrows)
      ba_vector <- numeric(nrows)

      survival_thin_ratio <- 1 # so that thin ratio isn't applied before thin age

      row <- 0

      # loop through all ages and calculate new values for n, ht, ba
      for (i in age1:age2) {
        if (i > age1) {
          n <- ((n^survival_params[1] + survival_params[2] * survival_thin_ratio *
                   (site_index/100)^survival_params[3] *
                   (i^survival_params[4] - (i-1)^survival_params[4]))^(1/survival_params[1]))
          ht <- predict_ht(age1 = i-1, age2 = i, ht = ht) *
            thin_response(age = i, thin_age = thin_age, thin_ratio = thin_ratio,
                          r = H_thin_params[1], k = H_thin_params[2]) * # additive assumption
            fert_response(age = i, fert_age = fert_age, N = kgNperha, P = P,
                          delta0 = H_fert_params[1], delta1 = H_fert_params[2],
                          beta = H_fert_params[3], gamma = H_fert_params[4])
          ba <- predict_ba(age1 = i-1, age2 = i, ba = ba) *
            thin_response(age = i, thin_age = thin_age, thin_ratio = thin_ratio,
                          r = BA_thin_params[1], k = BA_thin_params[2]) * # additive assumption
            fert_response(age = i, fert_age = fert_age, N = kgNperha, P = P,
                          delta0 = BA_fert_params[1], delta1 = BA_fert_params[2],
                          beta = BA_fert_params[3], gamma = 1)
        }
        row <- row + 1

        age_vector[row] <- i
        n_vector[row] <- n
        ht_vector[row] <- ht
        ba_vector[row] <- ba

        # Add extra row at thinning age to reduce basal area and surviving trees
        if (i == thin_age && isTRUE(valid_thin)) {
          ba <- ba * thin_ratio
          n <- n * thin_ratio # temporary fix, exploring models

          survival_thin_ratio <- thin_ratio # only apply after thin age

          row <- row + 1
          age_vector[row] <- i
          n_vector[row] <- n
          ht_vector[row] <- ht
          ba_vector[row] <- ba
        }
      }
      if(units == "imperial") {
        n_vector <- n_vector * acre_to_hectare
        ht_vector <- ht_vector / ft_to_m
        ba_vector <- ba_vector * m2_ha_to_ft2_ac
      }

      data.frame(
        runID = r,
        age = age_vector,
        n = round(n_vector, 0), # round to nearest whole tree
        ht = ht_vector,
        ba = ba_vector
      )

    })} # End loop

  out <- do.call(rbind, combined_data)

  return(out)
}
