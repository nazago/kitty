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

        ///////////ILP solver///////////
        
      
