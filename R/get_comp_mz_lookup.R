get_comp_mz_lookup <- function(.data, comp_name, rt_val, ppm){
  d = list()
  gln_mass = get_comp_mass(comp_name, mass_dict)
  c = get_element_count(comp_name)[['C']]
  n = get_element_count(comp_name)[['N']]
  comp_prefix <- .data %>%
    filter(formula == comp_name & rt == rt_val) %>%
    pull(prefix)
  ###two loops -
  ##one with compunds with just C and no N
  #one for compounds with C + N
  if(length(c) > 0)
  {
    #still calculate the mass of the compound
    #even if there are no Nitrogens and just look at all of the possible
    #carbon isotopes
    if(length(n) == 0)
    {
      #calculate all the possible carbon isotopes
      #and
      for(i in 0:c)
      {
        #number of Carbons and numbere of Nitrogens
        no_of_c = i
        no_of_n = 0
        #name the compound according to the 13C and 15N molecules
        compound = paste0(comp_prefix, '_', no_of_n, 'N', no_of_c, 'C')
        weight_compound = gln_mass + (i * w_c)
        #print(weight_compound)
        #print(" is weight_compound ")
        deviance <- ppm*0.000001*weight_compound
        #set upper and lower bounds on the masses
        #of the C and N isotopes
        upper_bound = weight_compound + deviance
        lower_bound = weight_compound - deviance
        #store, for each compound, the upper and lower bound masses
        #in a dictionary

        #if n is null, set to 0
        n = ifelse(is.null(n), 0, n)
        isotopeNumbers = n + c
        d[[compound]] = list('lb' = lower_bound,
                             'wc' = weight_compound,
                             'ub' = upper_bound,
                             'carbon' = no_of_c,
                             'nitrogen'= 0,
                             'isotope_numbers' =  isotopeNumbers)
      }
    }
    #if there are nitrogrens, no go ahead and calculate all masses
    #of all the possible C and N isotopes
    if(length(n) > 0)
    {
      #all possible combinations of 13C incorporation
      for(i in 0:c)
      {
        #possible combinations of 15N incorporation
        for(j in 0:n)
        {
          #number of Carbons and numbere of Nitrogens
          no_of_c = i
          no_of_n = j
          #name the compound according to the 13C and 15N molecules
          compound = paste0(comp_prefix, '_', no_of_n, 'N', no_of_c, 'C')
          weight_compound = gln_mass + (i * w_c) + (j * w_n)
          #print(weight_compound)
          #print("is the weight of compound")
          deviance <- ppm*0.000001*weight_compound
          #set upper and lower bounds on the masses
          #of the C and N isotopes
          upper_bound = weight_compound + deviance
          lower_bound = weight_compound - deviance
          #store, for each compound, the upper and lower bound masses
          #in a dictionary
          #if n is null, set to 0
          n = ifelse(is.null(n), 0, n)
          isotopeNumbers = n + c
          d[[compound]] = list('lb' = lower_bound,
                               'wc' = weight_compound,
                               'ub' = upper_bound,
                               'carbon' = no_of_c,
                               'nitrogen'= no_of_n,
                               'isotope_numbers' =  isotopeNumbers)

        }
      }
    }
  }
  return(d)
}
