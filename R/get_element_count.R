
#'get_element_count



#### INNITIALIZATION ####

w_c = 1.003355

w_n = 0.997035
ppm = 20
rt_tolerance = 0.1
mass_dict = list('N' = 14.003074,
                 'C' = 12.000000,
                 'O' = 15.994915,
                 'H' = 1.0078250,
                 'S' = 31.972072,
                 'P' = 30.973763,
                 'Cl' = 34.968853,
                 'Fe' = 55.934939,
                 'Cu' = 62.929599,
                 'Ca' = 39.962591,
                 'K' = 38.963708,
                 'Br' = 78.91833,
                 'F' = 18.998405 )

get_element_count <- function(chemical_formula){
  chemical_formula %>%
    makeup() %>%
    as.list() %>%
    return()
}
