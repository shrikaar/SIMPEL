get_comp_mass <- function(chemical_formula, mass_dict){
  #mass = mass_dict[['H']]
  mass = 1.007276
  element_list = get_element_count(chemical_formula)
  for(ele in names(element_list)){
    mass = mass + (element_list[[ele]]*mass_dict[[ele]])
  }
  return(mass)
}
