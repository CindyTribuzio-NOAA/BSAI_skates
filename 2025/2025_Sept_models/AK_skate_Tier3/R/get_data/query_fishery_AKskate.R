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
query_bsai_akskt <- function(new_year = 9999,
                             cas_sp = "USKT",
                             cas_sp_code = 703,
                             #fsh_target = "c",
                             cas_subarea = c("BS", "AI"),
                             obs_sp = 'ALASKA SKATE',
                             obs_sp_code = 88,
                             #twl_srvy = 98,
                             #srv_sp = 471,
                             area = 'bsai'){
  
  # get connected
  db = 'akfin'
  conn = afscdata::connect(db)  
  
# fishery data ----

## catch data ----
cat("\u231b", crayon::blue("working on catch query..."), "\n")
### query catch data and write raw data to folder ----
  #NOTE this still requires parsing of "skates, other" DO NOT USE in SS3 yet
#SimDesign::quiet(afscdata::q_catch(year = new_year,
#                                   species = cas_sp,
#                                   area = cas_subarea,
#                                   db = conn))

## fishery length data ----
cat("\u231b", crayon::blue("working on fishery length query..."), "\n")
dplyr::tbl(conn, dplyr::sql('norpac.debriefed_haul_mv')) %>% 
  dplyr::inner_join(dplyr::tbl(conn, dplyr::sql('norpac.debriefed_spcomp_mv')) %>% 
                      dplyr::filter(SPECIES == obs_sp_code),
                    by = c('HAUL_JOIN')) %>% 
  dplyr::inner_join(dplyr::tbl(conn, dplyr::sql('norpac.debriefed_length_mv')) %>% 
                      dplyr::filter(SPECIES == obs_sp_code),
                    by = c('HAUL_JOIN')) %>% 
  dplyr::rename_all(tolower) %>% 
  dplyr::select(gear,
                haul_join,
                numb = extrapolated_number,
                cruise = cruise.x,
                permit = permit.y,
                haul = haul.x,
                weight = extrapolated_weight,
                sex = sex.y,
                length,
                freq = frequency,
                lon = londd_end.x,
                lat = latdd_end.x,
                hday = haul_date.x,
                area = nmfs_area.x) %>% 
  dplyr::filter(area < 600) -> akfin_len 

dplyr::collect(akfin_len) %>% 
  dplyr::mutate(haul_join = paste0('H', haul_join),
                weight = weight / 1000,
                year = lubridate::year(hday),
                month = lubridate::month(hday),
                season = dplyr::case_when(month <= 2 ~ 1,
                                          month %in% c(3, 4) ~ 2,
                                          month %in% c(5, 6, 7, 8) ~ 3,
                                          month %in% c(9, 10) ~ 4,
                                          month >= 11 ~ 5),
                quarter = dplyr::case_when(month <= 3 ~ 1,
                                           month %in% c(4, 5, 6) ~ 2,
                                           month %in% c(7, 8, 9) ~ 3,
                                           month >= 10 ~ 4),
                trimester = dplyr::case_when(month <= 4 ~ 1,
                                             month %in% c(5, 6, 7, 8) ~ 2,
                                             month >= 9 ~ 3),
                gear = dplyr::case_when(gear %in% c(1, 2, 3, 4) ~ 'trawl',
                                        gear == 6 ~ 'pot',
                                        gear %in% c(5, 7, 9, 10, 11, 68, 8) ~ 'longline')) %>%  
  dplyr::filter(year <= new_year) %>% 
  vroom::vroom_write(., here::here(new_year, '2025_Sept_models', 'AK_skate_Tier3', 'data', 'raw', 'fish_lfreq.csv'), delim = ",")

capture.output(dplyr::show_query(akfin_len), 
               file = here::here(new_year, '2025_Sept_models', 'AK_skate_Tier3', "data", "sql", "fsh_lfreq_sql.txt"))
# print message when done
cat(crayon::green$bold("\u2713"), crayon::blue("fishery length query"), crayon::green$underline$bold$italic("DONE"), "\n")

## fishery observer catch data ----
cat("\u231b", crayon::blue("working on fishery observer catch query..."), "\n")
dplyr::tbl(conn, dplyr::sql('norpac.debriefed_haul_mv')) %>% 
  dplyr::inner_join(dplyr::tbl(conn, dplyr::sql('norpac.debriefed_spcomp_mv')) %>% 
                      dplyr::filter(SPECIES == obs_sp_code),
                    by = c('HAUL_JOIN')) %>% 
  dplyr::rename_all(tolower) %>% 
  dplyr::select(gear = gear_type,
                weight_kg = extrapolated_weight,
                lon = londd_start,
                lat = latdd_start,
                hday = haul_date.x,
                area = nmfs_area) %>% 
  dplyr::filter(area < 600) -> akfin_obsc

dplyr::collect(akfin_obsc) %>% 
  dplyr::mutate(year = lubridate::year(hday),
                gear = dplyr::case_when(gear %in% c(1, 2, 3, 4) ~ 'trawl',
                                        gear == 6 ~ 'pot',
                                        gear %in% c(5, 7, 9, 10, 11, 68, 8) ~ 'longline')) %>%  
  dplyr::filter(year <= new_year) %>% 
  vroom::vroom_write(., here::here(new_year, '2025_Sept_models', 'AK_skate_Tier3', 'data', 'raw', 'fish_obs_catch.csv'), delim = ",")

capture.output(dplyr::show_query(akfin_obsc), 
               file = here::here(new_year, '2025_Sept_models', 'AK_skate_Tier3', "data", "sql", "fsh_obs_catch_sql.txt"))
# print message when done
cat(crayon::green$bold("\u2713"), crayon::blue("fishery observer catch query"), crayon::green$underline$bold$italic("DONE"), "\n")

## specs ----
SimDesign::quiet(afscdata::q_specs(year = new_year,
                                   species = cas_sp,
                                   area = area,
                                   db = conn,
                                   save = TRUE))
cat(crayon::green$bold("\u2713"), crayon::blue("specs query"), crayon::green$underline$bold$italic("DONE"), "\n")

}
