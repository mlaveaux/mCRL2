// Author(s): Maurice Laveaux and Menno Bartels
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL_PBES_PBES_QUOTIENT_H
#define MCRL_PBES_PBES_QUOTIENT_H

#include "mcrl2/data/data_expression.h"
#include "mcrl2/pbes/pbes_expression.h"
#include "mcrl2/pbes/pbes_symmetry.h"
#include "mcrl2/utilities/indexed_set.h"

#include <boost/asio.hpp>
#include <boost/asio/buffer.hpp>
#include <boost/asio/read_until.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/process.hpp>
#include <boost/process/search_path.hpp>

namespace mcrl2::pbes_system {

class pbes_quotient
{
public:
    pbes_quotient(const detail::permutation& pi, const pbes& pbes, const std::string& gap_path)
    {
        if (!gap_path.empty())
        {
            gap_process = boost::process::child(gap_path,
                boost::process::args({"-E", "-q"}), 
                boost::process::std_in < input_stream,  
                boost::process::std_out > output_stream);
        
            if (pi.is_identity())
            {
                // Empty permutation, return immediately.
                std::string gap_input = "grp := Group(());\n";
                mCRL2log(log::debug) << "Setting symmetry group in GAP: " << gap_input;        
                input_stream << gap_input;
                input_stream.flush();

                std::string line;
                std::getline(output_stream, line);
                mCRL2log(log::debug) << "Received from GAP: " << line << std::endl;
            }
            else
            {        
                // Set the group in gap
                std::stringstream gap_input;
                gap_input << "grp := Group([";

                // Convert permutation to cycle notation
                int num_variables = pbes.initial_state().parameters().size();
                std::vector<bool> visited(num_variables, false);

                for (size_t i = 0; i < num_variables; ++i) {
                    if (!visited[i] && pi[i] != i) {
                        gap_input << "(";
                        
                        size_t current = i;
                        bool first_element = true;
                        do {
                            if (!first_element) {
                                gap_input << ",";
                            }
                            gap_input << (current + 1); // GAP uses 1-based indexing
                            visited[current] = true;
                            current = pi[current];
                            first_element = false;
                        } while (current != i);
                        
                        gap_input << ")";
                    }
                }

                gap_input << "]);\n";

                // Write to GAP process
                mCRL2log(log::debug) << "Setting symmetry group in GAP: " << gap_input.str();        
                input_stream << gap_input.str();
                input_stream.flush();

                std::string line;
                while (std::getline(output_stream, line)) {
                    if (line.find(']') != std::string::npos) {
                        break;
                    }
                }
                mCRL2log(log::debug) << "Received from GAP: " << line << std::endl;
            }
        }

    }

    /// Apply the quotienting to a propositional variable instantiation
    propositional_variable_instantiation apply(const propositional_variable_instantiation& pvi)
    {
        if (gap_process.running() == false)
        {
            return pvi;
        }

        mCRL2log(log::debug) << "Applying quotient to PVI: " << pvi << std::endl;

        m_temp_values.clear();
        for (const data::data_expression& param : pvi.parameters())
        {
            const auto& [index, inserted] = m_values.insert(param);
            m_temp_values.emplace_back(index);
        }

        std::stringstream gap_input;
        gap_input << "Minimum(List(Elements(grp), g -> Permuted([";
        for (size_t i = 0; i < m_temp_values.size(); ++i)
        {
            if (i > 0)
            {
                gap_input << ",";
            }
            gap_input << (m_temp_values[i] + 1); // GAP uses 1-based indexing
        }
        gap_input << "], g)));\n";

        mCRL2log(log::debug) << "Computing minimum using GAP: " << gap_input.str();
        input_stream << gap_input.str();
        input_stream.flush();

        // Read the result from GAP
        std::string line;
        std::string result_line;
        while (std::getline(output_stream, line)) {
            result_line += line;
            if (line.find(']') != std::string::npos) {
                break;
            }
        }
        line = result_line;

        mCRL2log(log::debug) << "Received from GAP: " << line << std::endl;

        // Parse the result
        new_params.clear();
        size_t pos = line.find('[');
        size_t end_pos = line.find(']');
        if (pos != std::string::npos && end_pos != std::string::npos && end_pos > pos)
        {
            std::string params_str = line.substr(pos + 1, end_pos - pos - 1);
            std::stringstream params_stream(params_str);
            std::string index_str;
            while (std::getline(params_stream, index_str, ','))
            {
                int index = std::stoi(index_str) - 1; // Convert back to 0-based indexing
                new_params.emplace_back(m_values.at(index));
            }
        }

        auto result = propositional_variable_instantiation(pvi.name(), data::data_expression_list(new_params.begin(), new_params.end()));        
        mCRL2log(log::debug) << "Resulting PVI: " << result << std::endl;
        return result;
    }

private:
    utilities::indexed_set<data::data_expression> m_values;

    std::vector<int> m_temp_values;
    std::vector<data::data_expression> new_params;

    boost::process::child gap_process;

    boost::process::ipstream output_stream;
    boost::process::opstream input_stream;
};

} // namespace mcrl2::pbes_system

#endif // MCRL_PBES_PBES_QUOTIENT_H