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
        
        for (uint8_t d = 0; d < tt.num_vars(); d++) {
        
         	  auto neg = cofactor0(tt, d);
		  auto pos = cofactor1(tt, d);
	
            if ( implies(pos, neg ) ) {
                flip_inplace(ttstar, d);
                inverted.at(d) = true; 
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


	auto ON_set = isop(ttstar);
        auto OFF_set = isop(unary_not(ttstar));
        
        
        set_add_rowmode(lp, TRUE); 

        
        
        for (cube cube : ON_set) { 

            for (uint8_t p = 0; p < tt.num_vars(); p++) { 

                auto no_literal_cube = cube;
                no_literal_cube.remove_literal(p);
		bool term = cube.num_literals() != no_literal_cube.num_literals();

                if (term) {
              
                    colno[p] = p + 1;
                    row[p] = 1;
                }
                else {
                    colno[p] = p + 1;
                    row[p] = 0;
                }
            }
            colno[tt.num_vars()] = tt.num_vars() + 1; 
            row[tt.num_vars()] = -1;            //
            add_constraintex(lp, tt.num_vars() + 1, row, colno, GE, 0);
        }


 
        for (auto& cube : OFF_set) { 

            for (uint8_t w = 0; w < tt.num_vars(); w++) { 

                auto no_literal_cube = cube;
                no_literal_cube.remove_literal(w);
		bool term = cube.num_literals() == no_literal_cube.num_literals();

                if (term) {
          
                    colno[w] = w + 1;
                    row[w] = 1;
                }
                else {
     
                    colno[w] = w + 1;
                    row[w] = 0;
                }
            }
            colno[tt.num_vars()] = tt.num_vars() + 1; 
            row[tt.num_vars()] = -1;            
            add_constraintex(lp, tt.num_vars() + 1, row, colno, LE, -1);
        }


        set_add_rowmode(lp, FALSE); 
     	
     	
     	for (uint8_t i = 0; i < tt.num_vars() + 1; i++) {
            
                colno[i] = i + 1;
                row[i] = 1;
            
            add_constraintex(lp, 1, row, colno, GE, 0);
        }
        
        for (uint8_t h = 0; h < tt.num_vars() + 1; h++) {
            colno[h] = h + 1;
            row[h] = 1;
        }
          
        set_obj_fnex(lp, tt.num_vars() + 1, row, colno);

        set_verbose(lp, IMPORTANT);
        set_minim(lp);
        int result = solve(lp);
        get_variables(lp, row);

     
        for(uint8_t u = 0; u < tt.num_vars(); u++) {
        	if (inverted[u]){
        		linear_form[u] = -row[u];
        		row[tt.num_vars()] = row[tt.num_vars()] -row[u];
        	}else{
        		linear_form[u] = row[u];
        		
        	}
            
            
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
