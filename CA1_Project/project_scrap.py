# Create a list of lists of all names that are not a part of any spine scaling list
no_spine_name_list2 = []
for comp in distList2[0]:
    if comp not in chain(list(chain(*no_spine_name_list)),
			 list(chain(*basal_name_list)),
			 list(chain(*med_spine_rad_name_list)),
			 list(chain(*med_spine_LM_name_list)),
			 list(chain(*max_spine_rad_name_list)),
			 list(chain(*thin_rad_name_list)),
			 list(chain(*thin_LM_name_list))):
        no_spine_name_list2.append(comp)


# Work in progress code for distance dependent conductance
#if (dist_dep == True):
#	    for comp in moose.wildcardFind(cell.path + '/' + '#[TYPE=' + comp_type + ']'):
#	        distance = np.sqrt(comp.x*comp.x + comp.y*comp.y + comp.z*comp.z)
#                for chan_name, cond in condSet.items():
#		    for dist,cond in condSet.items():
#		        if (dist[0] <= distance < dist[1]):                          
#                            conductance = cond
#		    addOneChan(library_name, chan_name, conductance, comp)   
            # Add the calcium pool to each compartment in the cell if it has been specified
#            if (CaPoolParams != None):
#	        add_calcium(library_name, cell, CaPoolParams, comp_type)
#	        for key in channelSet.keys():
#	            if ("Ca" in key):
#		        connect_cal2chan(channelSet[key].name, channelSet[key].chan_type, cell,
#                    		         CaPoolParams.caName, comp_type)
