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
        uint8_t variables_num = tt.num_vars();
        std::vector<int64_t> linear_form(variables_num + 1);
	std::vector<bool> complemented_vars(variables_num);
	TT TF_TruthTable = tt;
	
	

	for (uint8_t i = 0u; i < variables_num; i++)
	{
	uint8_t flag_n, flag_p =0;
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

        //ILP solver
        lprec* lp;
        int *colno = NULL;
        REAL *row = NULL;
        lp = make_lp(0, variables_num + 1);
        for (int i = 1; i <= (int) variables_num + 1; i++) { 
            //Mark all variables as integer numbers
            set_int(lp, i, TRUE);
        }
        colno = (int *) malloc((variables_num + 1) * sizeof(*colno));
        row = (REAL *) malloc((variables_num + 1) * sizeof(*row));

        set_add_rowmode(lp, TRUE); //The constraints are added row by row

        //onset
        std::vector<cube> onset_iso = isop(TF_TruthTable);
        std::vector<cube> offset_iso = isop(unary_not(TF_TruthTable));
        for (auto& c_cube : onset_iso) { //for every cube

            for (uint8_t var_index = 0; var_index < variables_num; var_index++) { //for every literal

                auto c_cube_without_literal = c_cube;
                c_cube_without_literal.remove_literal(var_index);

                //Add the constraint
                if (c_cube.num_literals() != c_cube_without_literal.num_literals()) {
                    //the cube contains the literal
                    colno[var_index] = var_index + 1;
                    row[var_index] = 1;
                }
                else {
        
                    colno[var_index] = var_index + 1;
                    row[var_index] = 0;
                }
            }
            colno[variables_num] = variables_num + 1; 
            row[variables_num] = -1;            
            add_constraintex(lp, variables_num + 1, row, colno, GE, 0);
        }

        //offset
        
        for (auto& c_cube : offset_iso) { 

            for (uint8_t var_index = 0; var_index < variables_num; var_index++) { //for every literal

                auto c_cube_without_literal = c_cube;
                c_cube_without_literal.remove_literal(var_index);

                //ILP solver
                if (c_cube.num_literals() == c_cube_without_literal.num_literals()) {
                   
                    colno[var_index] = var_index + 1;
                    row[var_index] = 1;
                }
                else {
                    
                    colno[var_index] = var_index + 1;
                    row[var_index] = 0;
                }
            }
            colno[variables_num] = variables_num + 1; 
            row[variables_num] = -1;            //
            add_constraintex(lp, variables_num + 1, row, colno, LE, -1);
        }

      
        for (uint8_t i = 0; i < variables_num + 1; i++) {
            for (uint8_t j = 0; j < variables_num + 1; j++) {
                colno[j] = j + 1;
                row[j] = (i == j) ? 1 : 0;
            }
            add_constraintex(lp, variables_num + 1, row, colno, GE, 0);
        }

        set_add_rowmode(lp, FALSE); 

       
        for (uint8_t var_index = 0; var_index < variables_num + 1; var_index++) {
            colno[var_index] = var_index + 1;
            row[var_index] = 1;
        }
        set_obj_fnex(lp, variables_num + 1, row, colno);

       
        set_verbose(lp, IMPORTANT);
        set_minim(lp);
        int result = solve(lp);
        get_variables(lp, row);

      
        for(uint8_t p = 0; p < variables_num; p++) 
        {
        if(complemented_vars[p])
        {
        linear_form[p] = -row[p];
        row[variables_num] = row[variables_num] - row[p];
        }
        else
        {
        linear_form[p] = row[p];
        row[variables_num] = row[variables_num];
        }
        }
        linear_form[variables_num] = (row[variables_num]);

        //Free allocated memory
        if(row != NULL)
            free(row);
        if(colno != NULL)
            free(colno);
        if(lp != NULL) 
            delete_lp(lp);

        //Return the results
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
