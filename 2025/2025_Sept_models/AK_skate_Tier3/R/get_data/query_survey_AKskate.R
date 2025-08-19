#' Query data for BSAI AK skate
#' Adapted by C Tribuzio from
#' Query data used for goa pcod
#' developed in 2024 by p. hulson
#' 
#' @param new_year current assessment year
#' @param fsh_sp species label for NORPAC data (default = 'ALASKA SKATE')
#' @param fsh_sp_code species code for NORPAC data (default = 88)
#' @param fsh_target trip target code for catch data (default = "c")
#' @param fsh_subarea NPFMC subareas (default to BSAI subareas)
#' @param twl_srvy gap survey id number (default = 98 for EBS shelf survey)
#' @param srv_sp gap species code (default = 471)
#' @param area assessment area (default = 'ebs_shelf')
#' 
query_ebssurv_akskt <- function(new_year = 9999,
                           #fsh_sp = "ALASKA SKATE",
                           #fsh_sp_code = 88,
                           #fsh_target = "c",
                           #fsh_subarea = c("EBS", "AI"),
                           twl_srvy = 98,
                           srv_sp = 471,
                           area = 'ebs_shelf'){
  
  # get connected
  db = 'akfin'
  conn = afscdata::connect(db)  
  
  # survey data ----
  
  ## trawl survey index data ----
  cat("\u231b", crayon::blue("working on trawl survey index query..."), "\n")
  dplyr::tbl(conn, dplyr::sql('gap_products.akfin_biomass')) %>% 
    dplyr::rename_all(tolower) %>% 
    dplyr::select(year,
                  survey_definition_id,
                  area_id,
                  species_code,
                  biomass_mt,
                  biomass_var,
                  population_count,
                  population_var) %>% 
    dplyr::filter(year <= new_year,
                  survey_definition_id == twl_srvy,
                  species_code == srv_sp,
                  area_id %in% c(99901)) %>% 
    dplyr::select(year,
                  survey = survey_definition_id, 
                  strata = area_id, 
                  species_code,
                  biom = biomass_mt, 
                  biom_var = biomass_var,
                  num = population_count,
                  num_var = population_var) %>%  
    dplyr::mutate(area = case_when(strata == 99901 ~ 'ebs_shelf')) -> twl_q
  
  dplyr::collect(twl_q) %>% 
    vroom::vroom_write(., here::here(new_year, '2025_Sept_models', 'data', 'raw', 'twl_srvy_index.csv'), delim = ",")
  capture.output(dplyr::show_query(twl_q), 
                 file = here::here(new_year, '2025_Sept_models', "data", "sql", "twl_srvy_index_sql.txt"))
  # print message when done
  cat(crayon::green$bold("\u2713"), crayon::blue("trawl survey index query"), crayon::green$underline$bold$italic("DONE"), "\n")
  
  ## trawl survey length pop'n data ----
  cat("\u231b", crayon::blue("working on trawl survey length pop'n query..."), "\n")
  dplyr::tbl(conn, dplyr::sql('gap_products.akfin_sizecomp')) %>% 
    dplyr::rename_all(tolower) %>% 
    dplyr::select(year,
                  survey_definition_id,
                  area_id,
                  species_code,
                  length_mm,
                  sex,
                  population_count) %>% 
    dplyr::filter(year <= new_year,
                  survey_definition_id == twl_srvy,
                  species_code == srv_sp,
                  area_id %in% c(99901)) %>% 
    dplyr::select(year,
                  survey = survey_definition_id, 
                  strata = area_id, 
                  species_code,
                  length = length_mm, 
                  sex,
                  num = population_count) -> twl_q
  
  dplyr::collect(twl_q) %>% 
    tidytable::filter(length > 0) %>% 
    tidytable::mutate(length = length / 10) %>% 
    vroom::vroom_write(., here::here(new_year, '2025_Sept_models', 'data', 'raw', 'twl_srvy_lpop.csv'), delim = ",")
  capture.output(dplyr::show_query(twl_q), 
                 file = here::here(new_year, '2025_Sept_models', "data", "sql", "twl_srvy_lpop_sql.txt"))
  # print message when done
  cat(crayon::green$bold("\u2713"), crayon::blue("trawl survey length pop'n query"), crayon::green$underline$bold$italic("DONE"), "\n")
  
  ## trawl survey age pop'n data ----
  cat("\u231b", crayon::blue("working on trawl survey age pop'n query..."), "\n")
  dplyr::tbl(conn, dplyr::sql('gap_products.akfin_agecomp')) %>% 
    dplyr::rename_all(tolower) %>% 
    dplyr::select(year,
                  survey_definition_id,
                  area_id,
                  species_code,
                  age,
                  sex,
                  population_count,
                  length_mm_mean) %>% 
    dplyr::filter(year <= new_year,
                  survey_definition_id == twl_srvy,
                  species_code == srv_sp,
                  area_id %in% c(99901)) %>% 
    dplyr::select(year,
                  survey = survey_definition_id, 
                  strata = area_id, 
                  species_code,
                  age, 
                  sex,
                  num = population_count,
                  mean_len = length_mm_mean) -> twl_q
  
  dplyr::collect(twl_q) %>% 
    tidytable::filter(age > 0) %>% 
    tidytable::mutate(mean_len = mean_len / 10) %>% 
    vroom::vroom_write(., here::here(new_year, '2025_Sept_models', 'data', 'raw', 'twl_srvy_apop.csv'), delim = ",")
  capture.output(dplyr::show_query(twl_q), 
                 file = here::here(new_year, '2025_Sept_models', "data", "sql", "twl_srvy_apop_sql.txt"))
  # print message when done
  cat(crayon::green$bold("\u2713"), crayon::blue("trawl survey age pop'n query"), crayon::green$underline$bold$italic("DONE"), "\n")
  
  ## trawl survey age specimen data ----
  cat("\u231b", crayon::blue("working on trawl survey age specimen query..."), "\n")
  dplyr::tbl(conn, dplyr::sql('gap_products.akfin_haul')) %>% 
    dplyr::inner_join(dplyr::tbl(conn, dplyr::sql('gap_products.akfin_cruise')),
                      by = c('CRUISEJOIN')) %>% 
    dplyr::inner_join(dplyr::tbl(conn, dplyr::sql('gap_products.akfin_specimen')),
                      by = c('HAULJOIN')) %>% 
    dplyr::rename_all(tolower) %>% 
    dplyr::select(year,
                  survey_definition_id,
                  species_code,
                  stratum,
                  hauljoin,
                  latitude_dd_start,
                  latitude_dd_end,
                  longitude_dd_start,
                  longitude_dd_end,
                  sex,
                  length_mm,
                  weight_g,
                  age) %>% 
    dplyr::filter(survey_definition_id %in% twl_srvy,
                  species_code %in% srv_sp) %>% 
    dplyr::mutate(lat_mid = (latitude_dd_start + latitude_dd_end) / 2,
                  long_mid = (longitude_dd_start + longitude_dd_end) / 2) %>% 
    dplyr::select(year, 
                  survey = survey_definition_id,
                  species_code,
                  stratum,
                  hauljoin,
                  sex,
                  length = length_mm,
                  weight = weight_g,
                  age,
                  lat_mid,
                  long_mid) -> twl_q
  
  dplyr::collect(twl_q) %>% 
    tidytable::filter(age > 0) %>% 
    tidytable::mutate(length = length / 10,
                      weight = weight / 1000) %>% 
    vroom::vroom_write(., here::here(new_year, '2025_Sept_models', 'data', 'raw', 'twl_srvy_age.csv'), delim = ",")
  capture.output(dplyr::show_query(twl_q), 
                 file = here::here(new_year, '2025_Sept_models', "data", "sql", "twl_srvy_age_sql.txt"))
  # print message when done
  cat(crayon::green$bold("\u2713"), crayon::blue("trawl survey age specimen query"), crayon::green$underline$bold$italic("DONE"), "\n")
  
  ## trawl survey age specimen data for age bias ----
  #dplyr::tbl(conn, dplyr::sql('gap_products.akfin_haul')) %>% 
  #  dplyr::inner_join(dplyr::tbl(conn, dplyr::sql('gap_products.akfin_cruise')),
  #                    by = c('CRUISEJOIN')) %>% 
  #  dplyr::inner_join(dplyr::tbl(conn, dplyr::sql('gap_products.akfin_specimen')),
  #                    by = c('HAULJOIN')) %>% 
  #  dplyr::rename_all(tolower) %>% 
  #  dplyr::select(survey_definition_id,
  #                species_code,
  #                cruise,
  #                vessel_id,
  #                haul,
  #                specimen_id,
  #                sex,
  #                length_mm,
  #                weight_g,
  #                age) %>% 
  #  dplyr::filter(survey_definition_id %in% c(47, 98),
  #                species_code %in% srv_sp,
  #                cruise %in% 200301,
  #                vessel_id %in% c(88, 89, 147)) %>% 
  #  dplyr::rename(survey = survey_definition_id,
  #                specimen = specimen_id,
  #                vessel = vessel_id,
  #                length = length_mm,
  #                weight = weight_g,
  #                original_age = age) -> twl_q
  #
  #dplyr::collect(twl_q) %>% 
  #  tidytable::filter(original_age > 0) %>% 
  #  tidytable::mutate(length = length / 10,
  #                    weight = weight / 1000) %>% 
  #  vroom::vroom_write(., here::here(new_year, 'data', 'raw', 'twl_srvy_age_bias.csv'), delim = ",")
  #capture.output(dplyr::show_query(twl_q), 
  #               file = here::here(new_year, "data", "sql", "twl_srvy_age_bias_sql.txt"))
  ## print message when done
  #cat(crayon::green$bold("\u2713"), crayon::blue("trawl survey age specimen for age bias query"), crayon::green$underline$bold$italic("DONE"), "\n")
  
  ## trawl survey haul info query ----
  cat("\u231b", crayon::blue("working on trawl survey haul data query..."), "\n")
  dplyr::tbl(conn, dplyr::sql('gap_products.akfin_haul')) %>% 
    dplyr::inner_join(dplyr::tbl(conn, dplyr::sql('gap_products.akfin_cruise')),
                      by = c('CRUISEJOIN')) %>%
    dplyr::rename_all(tolower) %>% 
    dplyr::select(year,
                  survey_definition_id,
                  stratum,
                  hauljoin,
                  gear_temperature_c,
                  surface_temperature_c,
                  depth_m,
                  latitude_dd_end,
                  longitude_dd_end) %>% 
    dplyr::filter(survey_definition_id %in% twl_srvy) %>% 
    dplyr::select(year, 
                  survey = survey_definition_id,
                  stratum,
                  hauljoin,
                  temp = gear_temperature_c,
                  stemp = surface_temperature_c,
                  depth = depth_m,
                  lat = latitude_dd_end,
                  lon = longitude_dd_end) -> twl_q
  
  dplyr::collect(twl_q) %>% 
    vroom::vroom_write(., here::here(new_year, '2025_Sept_models', 'data', 'raw', 'twl_srvy_haul.csv'), delim = ",")
  capture.output(dplyr::show_query(twl_q), 
                 file = here::here(new_year, '2025_Sept_models', "data", "sql", "twl_srvy_haul_sql.txt"))
  # print message when done
  cat(crayon::green$bold("\u2713"), crayon::blue("trawl survey haul data query"), crayon::green$underline$bold$italic("DONE"), "\n")
  
  ## trawl survey catch query ----
  cat("\u231b", crayon::blue("working on trawl survey catch data query..."), "\n")
  dplyr::tbl(conn, dplyr::sql('gap_products.akfin_haul')) %>% 
    dplyr::inner_join(dplyr::tbl(conn, dplyr::sql('gap_products.akfin_cruise')),
                      by = c('CRUISEJOIN')) %>% 
    dplyr::inner_join(dplyr::tbl(conn, dplyr::sql('gap_products.akfin_catch')),
                      by = c('CRUISEJOIN')) %>%
    dplyr::rename_all(tolower) %>% 
    dplyr::select(year,
                  survey_definition_id,
                  stratum,
                  hauljoin.x,
                  latitude_dd_end,
                  longitude_dd_end,
                  count,
                  distance_fished_km,
                  net_width_m,
                  species_code) %>% 
    dplyr::filter(survey_definition_id %in% twl_srvy,
                  species_code %in% srv_sp) %>% 
    dplyr::mutate(cpue = count / (distance_fished_km * (net_width_m * 1000))) %>% 
    dplyr::select(year, 
                  survey = survey_definition_id,
                  stratum,
                  hauljoin = hauljoin.x,
                  lat = latitude_dd_end,
                  lon = longitude_dd_end,
                  cpue,
                  species_code) -> twl_q
  
  dplyr::collect(twl_q) %>% 
    vroom::vroom_write(., here::here(new_year, '2025_Sept_models', 'data', 'raw', 'twl_srvy_catch.csv'), delim = ",")
  capture.output(dplyr::show_query(twl_q), 
                 file = here::here(new_year, '2025_Sept_models', "data", "sql", "twl_srvy_catch_sql.txt"))
  # print message when done
  cat(crayon::green$bold("\u2713"), crayon::blue("trawl survey catch data query"), crayon::green$underline$bold$italic("DONE"), "\n")
  
  ## trawl length freq query ----
  cat("\u231b", crayon::blue("working on trawl survey length freq data query..."), "\n")
  dplyr::tbl(conn, dplyr::sql('gap_products.akfin_haul')) %>% 
    dplyr::inner_join(dplyr::tbl(conn, dplyr::sql('gap_products.akfin_cruise')),
                      by = c('CRUISEJOIN')) %>% 
    dplyr::inner_join(dplyr::tbl(conn, dplyr::sql('gap_products.akfin_length')),
                      by = c('HAULJOIN')) %>%
    dplyr::rename_all(tolower) %>% 
    dplyr::select(year,
                  survey_definition_id,
                  stratum,
                  hauljoin,
                  gear_temperature_c,
                  surface_temperature_c,
                  depth_m,
                  latitude_dd_end,
                  longitude_dd_end,
                  species_code,
                  sex,
                  length_mm,
                  frequency) %>% 
    dplyr::filter(survey_definition_id %in% twl_srvy,
                  species_code %in% srv_sp) %>% 
    dplyr::mutate(length = length_mm / 10) %>% 
    dplyr::select(year, 
                  survey = survey_definition_id,
                  stratum,
                  hauljoin,
                  temp = gear_temperature_c,
                  stemp = surface_temperature_c,
                  depth = depth_m,
                  lat = latitude_dd_end,
                  lon = longitude_dd_end,
                  species_code,
                  sex,
                  length,
                  frequency) -> twl_q
  
  dplyr::collect(twl_q) %>% 
    vroom::vroom_write(., here::here(new_year, '2025_Sept_models', 'data', 'raw', 'twl_srvy_lenfreq.csv'), delim = ",")
  capture.output(dplyr::show_query(twl_q), 
                 file = here::here(new_year, '2025_Sept_models', "data", "sql", "twl_srvy_lenfreq_sql.txt"))
  # print message when done
  cat(crayon::green$bold("\u2713"), crayon::blue("trawl survey length freq data query"), crayon::green$underline$bold$italic("DONE"), "\n")

  ## trawl survey strata data ----
  cat("\u231b", crayon::blue("working on trawl survey strata data query..."), "\n")
  # strata with area sizes
  dplyr::tbl(conn, dplyr::sql('gap_products.akfin_area')) %>% 
    dplyr::rename_all(tolower) %>% 
    dplyr::filter(survey_definition_id %in% twl_srvy,
                  area_type == 'STRATUM') %>% 
    dplyr::select(survey = survey_definition_id,
                  design_year,
                  stratum = area_id,
                  area = area_km2) %>% 
    dplyr::collect() -> st_area
  # subregion level with description (e.g., wgoa, etc)
  dplyr::tbl(conn, dplyr::sql('gap_products.akfin_area')) %>% 
    dplyr::rename_all(tolower) %>% 
    dplyr::filter(survey_definition_id %in% twl_srvy,
                  area_type == 'REGULATORY AREA') %>% 
    dplyr::select(area_id,
                  subarea_name = description,
                  design_year) %>% 
    dplyr::collect() -> subreg
  # strata within subregions
  dplyr::tbl(conn, dplyr::sql('gap_products.akfin_stratum_groups')) %>% 
    dplyr::rename_all(tolower) %>% 
    dplyr::filter(survey_definition_id %in% twl_srvy) %>% 
    dplyr::select(stratum, 
                  area_id) %>% 
    dplyr::collect() -> st_subreg
  # join all to get strata with area sizes and subregion ids
  st_area %>% 
    tidytable::left_join(st_subreg %>% 
                           tidytable::left_join(subreg) %>% 
                           tidytable::drop_na()) %>%
    tidytable::filter(design_year == max(design_year), .by = c(stratum)) %>% 
    vroom::vroom_write(., here::here(new_year, 'data', 'raw', 'twl_srvy_strata.csv'), delim = ",")
  # print message when done
  cat(crayon::green$bold("\u2713"), crayon::blue("trawl survey strata data query"), crayon::green$underline$bold$italic("DONE"), "\n")

  ## longline survey index data ----
  #cat("\u231b", crayon::blue("working on longline survey index query..."), "\n")
  #SimDesign::quiet(afscdata::q_lls_rpn(year = new_year,
  #                    species = srv_sp,
  #                    area = area,
  #                    by = 'geoarea',
  #                    use_historical = FALSE,
  #                    db = conn))
  # print message when done
  #cat(crayon::green$bold("\u2713"), crayon::blue("longline survey index query"), crayon::green$underline$bold$italic("DONE"), "\n")
  
  ## longline survey length data ----
  #cat("\u231b", crayon::blue("working on longline survey length rpn query..."), "\n")
  #SimDesign::quiet(afscdata::q_lls_rpn_length(year = new_year,
  #                           species = srv_sp,
  #                           area = area,
  #                           db = conn))
  # print message when done
  #cat(crayon::green$bold("\u2713"), crayon::blue("longline survey length rpn query"), crayon::green$underline$bold$italic("DONE"), "\n")

  ## iphc survey index data ----
  #cat("\u231b", crayon::blue("working on iphc survey index query..."), "\n")
  #dplyr::tbl(conn, dplyr::sql('afsc_host.fiss_rpn')) %>% 
  #  dplyr::rename_all(tolower) %>% 
  #  dplyr::filter(species %in% c('Pacific cod'),
  #                fmp_sub_area %in% c("CGOA", "EY/SE", "WGOA", "WY")) %>% 
  #  dplyr::select(year = survey_year, 
  #                fmp_sub_area, 
  #                species, 
  #                strata = rpn_strata, 
  #                strata_rpn, 
  #                boot_sd, 
  #                boot_bias,
  #                n_stations,
  #                n_pos_catch) -> iphc_q
  #
  #dplyr::collect(iphc_q) %>% 
  #  vroom::vroom_write(., here::here(new_year, 'data', 'raw', 'iphc_srvy_index.csv'), delim = ",")
  #capture.output(dplyr::show_query(iphc_q), 
  #               file = here::here(new_year, "data", "sql", "iphc_srvy_index_sql.txt"))
  # print message when done
  #cat(crayon::green$bold("\u2713"), crayon::blue("iphc survey index query"), crayon::green$underline$bold$italic("DONE"), "\n")

}
