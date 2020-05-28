#this function is going to evaluate all of the compounds in the metadata
#id'd by mass and retention time and determine which of their possible
#C and N isotopes are present in the data
get_comp_stage <- function(x, y, comp_lookup_table, r_time, rt_tolerance){
  #include all the formula info for the compound as well
  d = list()
  for(k in names(comp_lookup_table)){
    lb = comp_lookup_table[[k]]['lb']
    ub = comp_lookup_table[[k]]['ub']
    c = comp_lookup_table[[k]]['carbon']
    n = comp_lookup_table[[k]]['nitrogen']
    #for a given compound in the lookup table
    #determine if its characteristics fall within the lower and upper
    #bounds of the mass-rt features actually in the data
    if (lb <= x & x <= ub){
      if((r_time - rt_tolerance) <= y & y <= (r_time + rt_tolerance)){

        d = list('compound' = k,
                 'carbon' = c,
                 'nitrogen' = n)
        return(d)
      }
    }
  }
}
