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
        std::vector<bool> inverted(tt.num_vars(), false);    
        auto ttstar = tt;    
        
        for (uint8_t var = 0; var < tt.num_vars(); var++) {
        
         	  auto neg = cofactor0(tt, var);
		  auto pos = cofactor1(tt, var);
	
            if ( implies(pos, neg ) ) {
                flip_inplace(ttstar, var);
                inverted.at(var) = true; 
            }
            else if  ( !implies(neg, pos) ) {

                return false;
            }
        }

    
        lprec* lp;
        int *colno = NULL;
        REAL *row = NULL;
        lp = make_lp(0, tt.num_vars() + 1);
        
        
        for ( uint64_t i = 1; i < tt.num_vars() + 1; i++) { 
         
            set_int(lp, i, TRUE);
        }
        
        
        colno = (int *) malloc((tt.num_vars() + 1) * sizeof(*colno));
        row = (REAL *) malloc((tt.num_vars() + 1) * sizeof(*row));


	std::vector<cube> onset_isop = isop(ttstar);
        std::vector<cube> offset_isop = isop(unary_not(ttstar));
        
        
        set_add_rowmode(lp, TRUE); 

        
        
        for (cube cube : onset_isop) { 

            for (uint8_t var = 0; var < tt.num_vars(); var++) { 

                auto cube_without_literal = cube;
                cube_without_literal.remove_literal(var);


                if (cube.num_literals() != cube_without_literal.num_literals()) {
              
                    colno[var] = var + 1;
                    row[var] = 1;
                }
                else {
                    colno[var] = var + 1;
                    row[var] = 0;
                }
            }
            colno[tt.num_vars()] = tt.num_vars() + 1; 
            row[tt.num_vars()] = -1;            //
            add_constraintex(lp, tt.num_vars() + 1, row, colno, GE, 0);
        }


 
        for (auto& cube : offset_isop) { 

            for (uint8_t var = 0; var < tt.num_vars(); var++) { 

                auto cube_without_literal = cube;
                cube_without_literal.remove_literal(var);


                if (cube.num_literals() == cube_without_literal.num_literals()) {
          
                    colno[var] = var + 1;
                    row[var] = 1;
                }
                else {
     
                    colno[var] = var + 1;
                    row[var] = 0;
                }
            }
            colno[tt.num_vars()] = tt.num_vars() + 1; 
            row[tt.num_vars()] = -1;            
            add_constraintex(lp, tt.num_vars() + 1, row, colno, LE, -1);
        }


      

        set_add_rowmode(lp, FALSE); 
     
        for (uint8_t var = 0; var < tt.num_vars() + 1; var++) {
            colno[var] = var + 1;
            row[var] = 1;
        }
        
          for (uint8_t i = 0; i < tt.num_vars() + 1; i++) {
            for (uint8_t j = 0; j < tt.num_vars() + 1; j++) {
                colno[j] = j + 1;
                row[j] = (i == j) ? 1 : 0;
            }
            add_constraintex(lp, tt.num_vars() + 1, row, colno, GE, 0);
        }
        set_obj_fnex(lp, tt.num_vars() + 1, row, colno);

        set_verbose(lp, IMPORTANT);
        set_minim(lp);
        int result = solve(lp);
        get_variables(lp, row);

     
        for(uint8_t var = 0; var < tt.num_vars(); var++) {
            linear_form[var] = (inverted[var]) ? -row[var] : row[var];
            row[tt.num_vars()] = (inverted[var]) ? row[tt.num_vars()] - row[var] : row[tt.num_vars()];
        }
        linear_form[tt.num_vars()] = (row[tt.num_vars()]);

        
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

} /* namespace kitty */ 
