#pragma once

#include <iostream>
#include "isop.hpp"
#include "cube.hpp"
#include <vector>
#include <lpsolve/lp_lib.h> /* uncomment this line to include lp_solve */
#include "traits.hpp"

namespace kitty {


    template<typename TT, typename = std::enable_if_t<is_complete_truth_table<TT>::value>>
    bool is_threshold(const TT& tt, std::vector<int64_t>* plf = nullptr)
    {
      
        std::vector<int64_t> linear_form(tt.num_vars() + 1);
	std::vector<bool> complemented_vars(tt.num_vars());
	TT TF_TruthTable = tt;
	std::vector<cube> On_set_vector = isop(TF_TruthTable);
        std::vector<cube> Off_set_vector = isop(unary_not(TF_TruthTable));
	

	for (uint8_t i = 0u; i < tt.num_vars(); i++)
	{
	uint64_t flag_n, flag_p =0;
	TT f0 = cofactor0(TF_TruthTable, i);
	TT f1 = cofactor1(TF_TruthTable, i);


	if (implies(f1,f0))
	{
		flag_n = 1;
	}
	else if (implies(f0,f1))
	{
		flag_p = 1;
	} 
	else
	{
	return false;
	} 
	

	if (flag_n == 1)
	{
		flip_inplace( TF_TruthTable, i );
		complemented_vars.at(i) = true;
	}   
	if (flag_p == 1)
	{
		complemented_vars.at(i) = false;
	}  

	} 

 
        
        
        
        lprec* lp;
        int *colno = NULL;
        REAL *row = NULL;
        lp = make_lp(0, tt.num_vars() + 1);
        for (int w = 1; w <= (int) tt.num_vars() + 1; w++) 
        { 
 
            set_int(lp, w, TRUE);
        }
        colno = (int *) malloc((tt.num_vars() + 1) * sizeof(*colno));
        row = (REAL *) malloc((tt.num_vars() + 1) * sizeof(*row));
        if((colno == NULL) || (row == NULL))
       	 return false;
        set_add_rowmode(lp, TRUE); 


        for (auto& ccube : On_set_vector) 
        { 

            for (uint8_t k = 0; k < tt.num_vars(); k++) 
            { 
                auto ccube_without_literal = ccube;
                ccube_without_literal.remove_literal(k);

               
                 if (ccube.num_literals() != ccube_without_literal.num_literals())
                 {
                   
                    colno[k] = k + 1;
                    row[k] = 1;
                 }
                 else 
                 {
                    
                    colno[k] = k + 1;
                    row[k] = 0;
                 }
            }
            colno[tt.num_vars()] = tt.num_vars() + 1; 
            row[tt.num_vars()] = -1;            
            add_constraintex(lp, tt.num_vars() + 1, row, colno, GE, 0);
        }

       
        
        for (auto& ccube : Off_set_vector) 
        { 

            for (uint64_t m = 0; m < tt.num_vars(); m++) 
            {

                auto ccube_without_literal = ccube;
                ccube_without_literal.remove_literal(m);

                
                if (ccube.num_literals() == ccube_without_literal.num_literals()) 
                {
                    
                    colno[m] = m + 1;
                    row[m] = 1;
                }
                else 
                {
                    
                    colno[m] = m + 1;
                    row[m] = 0;
                }
            }
            colno[tt.num_vars()] = tt.num_vars() + 1; 
            row[tt.num_vars()] = -1;            
            add_constraintex(lp, tt.num_vars() + 1, row, colno, LE, -1);
        }

        
       

        set_add_rowmode(lp, FALSE); 

        
        for (uint64_t l = 0; l < tt.num_vars() + 1; l++) 
        {
            colno[l] = l + 1;
            row[l] = 1;
        }
        set_obj_fnex(lp, tt.num_vars() + 1, row, colno);

        
        set_verbose(lp, IMPORTANT);
        set_minim(lp);
        int result = solve(lp);
        get_variables(lp, row);

        
        for(uint64_t p = 0; p < tt.num_vars(); p++) 
        {
        if(complemented_vars[p])
        {
        linear_form[p] = -row[p];
        row[tt.num_vars()] = row[tt.num_vars()] - row[p];
        }
        else
        {
        linear_form[p] = row[p];
        }
        }
        linear_form[tt.num_vars()] = (row[tt.num_vars()]);

        
        if(row != NULL)
        {
            free(row);
            }
        if(colno != NULL)
        {
            free(colno);
            }
        if(lp != NULL)
        { 
            delete_lp(lp);
            }

        
        if(result == 0) {
            if(plf)
            {
                *plf = linear_form;
                }
            return true;
        }
        else {
            return false;
        }
    }
    }
  
    
