#pragma once
#include <iostream>
#include "isop.hpp"
#include "cube.hpp"
#include "constructors.hpp"
#include "dynamic_truth_table.hpp"
#include "static_truth_table.hpp"
#include <vector>
#include <lpsolve/lp_lib.h> /* uncomment this line to include lp_solve */
#include "traits.hpp"

namespace kitty 
{

    template<typename TT, typename = std::enable_if_t<is_complete_truth_table<TT>::value>>
    bool is_threshold(const TT& tt, std::vector<int64_t>* plf = nullptr)
    {
        auto variables_num = tt.num_vars();
        std::vector<int64_t> linear_form(variables_num + 1);
	std::vector<bool> complemented_vars(variables_num);
	TT TF_TruthTable = tt;
	std::vector<cube> On_set_vector = isop(TF_TruthTable);
        std::vector<cube> Off_set_vector = isop(unary_not(TF_TruthTable));
	

	for (auto i = 0u; i < variables_num; i++)
	{
	uint8_t flag_n, flag_p =0;
	auto f0 = cofactor0(TF_TruthTable, i);
	auto f1 = cofactor1(TF_TruthTable, i);

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

        ///////////ILP solver///////////
        
        
        
        lprec* lp;
        int *colno = NULL;
        REAL *row = NULL;
        lp = make_lp(0, variables_num + 1);
       
        colno = (int *) malloc((variables_num + 1) * sizeof(*colno));
        row = (REAL *) malloc((variables_num + 1) * sizeof(*row));
        if((colno == NULL) || (row == NULL))
        return false;
        set_add_rowmode(lp, TRUE); 

        
        
        for (auto& ccube : On_set_vector) { 

            for (uint8_t k = 0; k < variables_num; k++) { 
                auto ccube_without_literal = ccube;
                ccube_without_literal.remove_literal(k);

               
                if (ccube.num_literals() != ccube_without_literal.num_literals()) {
                   
                    colno[k] = k + 1;
                    row[k] = 1;
                }
                else {
                    
                    colno[k] = k + 1;
                    row[k] = 0;
                }
            }
            colno[variables_num] = variables_num + 1; 
            row[variables_num] = -1;            
            add_constraintex(lp, variables_num + 1, row, colno, GE, 0);
        }

       
        
        for (auto& ccube : Off_set_vector) { 

            for (uint8_t k = 0; k < variables_num; k++) {

                auto ccube_without_literal = ccube;
                ccube_without_literal.remove_literal(k);

                
                if (ccube.num_literals() == ccube_without_literal.num_literals()) {
                    
                    colno[k] = k + 1;
                    row[k] = 1;
                }
                else {
                    
                    colno[k] = k + 1;
                    row[k] = 0;
                }
            }
            colno[variables_num] = variables_num + 1; 
            row[variables_num] = -1;            
            add_constraintex(lp, variables_num + 1, row, colno, LE, -1);
        }

        
       

        set_add_rowmode(lp, FALSE); 

        
        for (uint8_t k = 0; k < variables_num + 1; k++) {
            colno[k] = k + 1;
            row[k] = 1;
        }
        set_obj_fnex(lp, variables_num + 1, row, colno);

        
        set_verbose(lp, IMPORTANT);
        set_minim(lp);
        int result = solve(lp);
        get_variables(lp, row);

        
        for(uint8_t k = 0; k < variables_num; k++) 
        {
        if(complemented_vars[k])
        {
        linear_form[k] = -row[k];
        row[variables_num] = row[variables_num] - row[k];
        }
        else
        {
        linear_form[k] = row[k];
        row[variables_num] = row[variables_num];
        }
        }
        linear_form[variables_num] = (row[variables_num]);

        
        if(row != NULL)
            free(row);
        if(colno != NULL)
            free(colno);
        if(lp != NULL) 
            delete_lp(lp);

        
        if(result == 0) {
            if(plf)
                *plf = linear_form;
            return true;
        }
        else {
            return false;
        }
    }
    }
    
