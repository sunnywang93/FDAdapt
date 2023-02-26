
# Import dataset
powerconsumption_full <- read.table(unz("./data-raw/household_power_consumption.zip",
                                         "household_power_consumption.txt"),
                                    header = TRUE, sep = ";")


# Sort curves by Time, starting from midnight (00:00h)
powerconsumption_df <- powerconsumption_full |>
  dplyr::mutate(time_stamp = lubridate::hms(Time),
                hour = lubridate::hour(time_stamp),
                minute = lubridate::minute(time_stamp),
                second = lubridate::second(time_stamp)) |>
  dplyr::arrange(hour, minute, second) |>
  dplyr::select(Date, Time, Voltage) |>
  tibble::rowid_to_column("id")

#Convert to matrix
powerconsumption <- powerconsumption_df |>
  tidyr::pivot_wider(id_cols = -id,
                     names_from = Time,
                     values_from = Voltage) |>
  dplyr::select(-Date) |>
  as.matrix() |>
  unname()

#Convert from characters to numeric and delete rows with NAs
powerconsumption <- apply(powerconsumption, 2, as.numeric)
powerconsumption <- powerconsumption[complete.cases(powerconsumption),] |>
  tibble::as_tibble()

usethis::use_data(powerconsumption, overwrite = TRUE)





