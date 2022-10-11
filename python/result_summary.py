"""
Combine analysis results

Authors:
    Sergei L Kosakovsky Pond (spond@temple.edu)
    Alexander G Lucaci (alexander.lucaci@temple.edu)

Version:
    v0.0.1 (2021-01-17)
    v0.0.2 (2022-09-19)

"""

# #######################################################################
# Imports
# #######################################################################

import argparse
import csv
import random
import os
import json
import sys
import re
import math
import numpy
import progressbar
from tqdm import tqdm
from collections import Counter

# #######################################################################
# Declares
# #######################################################################
#bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
#bar = progressbar.ProgressBar()

random.seed ()
by_file = {}
timer = 0
count = 0
tags = {}
nucs = set (['A','C','G','T'])

# #######################################################################
# Command line arguments
# #######################################################################
arguments = argparse.ArgumentParser(description='Combine alignments into a single file, adding a reference sequence as well')
arguments.add_argument('-i', '--input',  help = 'Directories ', required = True, type = str, nargs = '*')
arguments.add_argument('-o', '--output',  help = 'Output file', required = True, type = str)
arguments.add_argument('-s', '--simulations',  help = 'Directories ', action = 'store_true' )
arguments.add_argument('-l', '--latex',  help = 'Make LATeX output ', action = 'store_true' )
settings = arguments.parse_args()

# #######################################################################
# Helper functions
# #######################################################################

def get_omega3 (fit):
    omegas = fit["fits"]["Unconstrained model"]["Rate Distributions"]["Test"]
    return omegas[str (len (omegas) - 1)]

def get_omega2 (fit):
    omegas = fit["fits"]["Unconstrained model"]["Rate Distributions"]["Test"]
    return omegas[str (len (omegas) - 2)]

def get_omega1 (fit):
    omegas = fit["fits"]["Unconstrained model"]["Rate Distributions"]["Test"]
    return omegas[str (0)]

def get_srv (fit):
    try:
        srv = fit["fits"]["Unconstrained model"]["Rate Distributions"]["Synonymous site-to-site rates"]
        s  = 0
        s2 = 0
        for k, v in srv.items():
            s += v["rate"] * v["proportion"] 
            s2 += v["rate"]*v["rate"] * v["proportion"] 
        return math.sqrt ((s2-s*s)/s)
    except:
        return 0

# #######################################################################
# Newick parser
# #######################################################################

def newick_parser(nwk_str, bootstrap_values, track_tags):
    clade_stack = []
    automaton_state = 0
    current_node_name = ""
    current_node_attribute = ""
    current_node_annotation = ""
    quote_delimiter = None
    name_quotes = {
      "'": 1,
      '"': 1
    }
    
    def add_new_tree_level():
      new_level = {
        "name": None
      };
      the_parent = clade_stack[len(clade_stack) - 1]
      if (not "children" in the_parent):
        the_parent["children"] = [];
      
      clade_stack.append (new_level);
      the_parent["children"].append(clade_stack[len(clade_stack) - 1]);
      clade_stack[len(clade_stack)-1]["original_child_order"] = len(the_parent["children"])
    

    def finish_node_definition():
      nonlocal current_node_name
      nonlocal current_node_annotation
      nonlocal current_node_attribute
      
      this_node = clade_stack.pop()
      if (bootstrap_values and "children" in this_node):
        this_node["bootstrap_values"] = current_node_name
      else:
        this_node["name"] = current_node_name
      
      this_node["attribute"] = current_node_attribute
      this_node["annotation"] = current_node_annotation
      
      try:
      
          if not 'children' in this_node:
            node_tag = "background"
            for k, v in tags.items():
                if this_node["name"].find (k) >= 0:
                    node_tag = v
                    break
          else:
            '''
            counts = {}
            node_tag = ""
            for n in this_node['children']:
                counts[n["tag"]] = 1 + (counts[n["tag"]] if n["tag"] in counts  else 0)
            if len (counts) == 1:
                node_tag = list (counts.keys())[0]
            '''
            node_tag = "test"
        
          this_node["tag"] = node_tag
      except Exception as e:
        print ("Exception ", e)
        
      if track_tags is not None:
        track_tags[this_node["name"]] = [this_node["tag"], 'children' in this_node]
       
      current_node_name = ""
      current_node_attribute = ""
      current_node_annotation = ""
    

    def generate_error(location):
      return {
        'json': None,
        'error':
          "Unexpected '" +
          nwk_str[location] +
          "' in '" +
          nwk_str[location - 20 : location + 1] +
          "[ERROR HERE]" +
          nwk_str[location + 1 : location + 20] +
          "'"
      }

    tree_json = {
      "name" : "root"
    }
    
    clade_stack.append(tree_json);

    space = re.compile("\s")

    for char_index in range (len(nwk_str)):
      try:
        current_char = nwk_str[char_index]
        if automaton_state == 0:
           #look for the first opening parenthesis
           if (current_char == "("):
              add_new_tree_level()
              automaton_state = 1
        elif automaton_state == 1 or automaton_state == 3:
            #case 1: // name
            #case 3: { // branch length
            #reading name
            if (current_char == ":"):
              automaton_state = 3;
            elif current_char == "," or current_char == ")":
              try:
                finish_node_definition()
                automaton_state = 1
                if (current_char == ","):
                  add_new_tree_level()
              except Exception as e:
                return generate_error(char_index)
              
            elif (current_char == "("):
              if len(current_node_name) > 0:
                return generate_error(char_index);
              else:
                add_new_tree_level()
              
            elif (current_char in name_quotes):
              if automaton_state == 1 and len(current_node_name) == 0 and len (current_node_attribute) == 0 and len (current_node_annotation) == 0:
                automaton_state = 2
                quote_delimiter = current_char
                continue
              return generate_error(char_index)
            else:
              if (current_char == "{"):
                if len (current_node_annotation):
                  return generate_error(char_index)
                else:
                  automaton_state = 4
              else:
                if (automaton_state == 3):
                  current_node_attribute += current_char;
                else:
                  if (space.search(current_char)):
                    continue;
                  if (current_char == ";"):
                    char_index = len(nwk_str)
                    break
                  current_node_name += current_char;
        elif automaton_state == 2: 
            # inside a quoted expression
            if (current_char == quote_delimiter):
              if (char_index < len (nwk_str - 1)):
                if (nwk_str[char_index + 1] == quote_delimiter):
                  char_index+=1
                  current_node_name += quote_delimiter;
                  continue;

              quote_delimiter = 0
              automaton_state = 1
              continue
            else:
              current_node_name += current_char;
        elif automaton_state == 4:
           ##inside a comment / attribute
            if (current_char == "}"):
              automaton_state = 3
            else:
              if (current_char == "{"):
                return generate_error(char_index);
              current_node_annotation += current_char;
      except Exception as e:
        return generate_error(char_index);

    if (len (clade_stack) != 1):
      return generate_error(len (nwk_str) - 1);

    if (len (current_node_name)):
        tree_json['name'] = current_node_name;

    return {
      'json': tree_json,
      'error': None
    }
#end method
    
def traverse_tree (node, parent, labels, counter, lengths, lookup, model):
    if node["name"] in labels:
        node["label"] = labels[node["name"]]
        if parent:
            diff = 0
            for i, c in enumerate (node["label"]):
                if c in nucs:
                    if parent["label"][i] != c and parent["label"][i] in nucs:
                        diff += 1
                        
            if diff > 0:
                counter[diff] += 1
                l = lookup[node["name"]][model]
                lengths[diff-1].append (l)
            
            
            counter[0] += 1
            
    else:
        node["label"] = parent["label"]
        counter[0] += 1
        
    if "children" in node:
        for c in node["children"]:
           traverse_tree (c, node, labels, counter, lengths, lookup, model) 
#end method 

# #######################################################################
# Main
# #######################################################################
i = 0
    
for basedir in settings.input:
     for root, dirs, files in os.walk(basedir):
        for each_file in tqdm(files):
            file_name, file_ext =  os.path.splitext (each_file)
            #print (file_ext, file = sys.stderr)
            if file_ext == '.json':
                #print(file_name)
                parts = file_name.split ('.')
                #print (parts)
                with open (os.path.join (root, each_file), "r") as fh:
                    try:
                        i += 1
                        #bar.update(i)
                        results = json.load (fh)
                        pv = results['test results']['p-value']
                        lrt = results['test results']['LRT']
                        logl = results["fits"]["Unconstrained model"]["Log Likelihood"]
                        aic = results["fits"]["Unconstrained model"]["AIC-c"]
                        omega3 = get_omega3
                        file_key = parts[0]
                        if settings.simulations:
                            for k in parts[1:]:
                                if k != "BUSTED":
                                    file_key += "_" + k
                                else:
                                    break
                        if not file_key in by_file:
                            by_file[file_key] = {}
                        model_key = parts[-1]
                        by_file[file_key][model_key] = {'p' : pv, 
                                                        'LR' : lrt, 
                                                        'logL' : logl, 
                                                        'AIC' : aic, 
                                                        'omega'  : get_omega3 (results),
                                                        'omega1' : get_omega1 (results),
                                                        'omega2' : get_omega2 (results),
                                                        'SRV' : get_srv (results)}
                                                        
                        if model_key == "MH":
                            tree = newick_parser (results["input"]["trees"]["0"],{},{})["json"]
                            counter = [0,0,0,0]
                            lengths = [[],[],[]]
                            
                             
                            for site, data in results["substitutions"]["0"].items():
                                traverse_tree (tree,  None, data, counter, lengths, results['branch attributes']['0'], 'unconstrained')

                            
                            lengths = ([numpy.mean (l) if len (l) else 0. for l in lengths])
                            by_file[file_key][model_key]['subs'] = counter
                            by_file[file_key][model_key]['N'] = results['input']['number of sequences']
                            by_file[file_key][model_key]['S'] = results['input']['number of sites']
                            by_file[file_key][model_key]['T'] = sum([v['unconstrained'] if 'unconstrained' in v else 0 for k,v in results['branch attributes']['0'].items()])
                            by_file[file_key][model_key]['delta'] = results["fits"]["Unconstrained model"]["Rate Distributions"]["rate at which 2 nucleotides are changed instantly within a single codon"]
                            by_file[file_key][model_key]['psi'] = results["fits"]["Unconstrained model"]["Rate Distributions"]["rate at which 3 nucleotides are changed instantly within a single codon"]
                            prior = by_file[file_key][model_key]['omega']
                            by_file[file_key][model_key]['L'] = lengths
                            if prior ['proportion'] == 1: 
                                prior = 0
                            else:
                                prior = prior ['proportion'] / (1 -  prior ['proportion'])
                            EB = 0
                            if prior > 0:   
                                for b, data in results['branch attributes']['0'].items():
                                    if "Posterior prob omega class by site" in data:
                                        for s in data["Posterior prob omega class by site"][2]:
                                            if s > 0 and s < 1 and s/(1-s) / prior >= 100:
                                                EB+=1
                            by_file[file_key][model_key]['EB'] = EB
                                    
                    except Exception as e:
                        print (e, each_file, file = sys.stderr)
                        raise
                        pass
                    
print("# Done processing")

# #######################################################################
# LaTex Output or normal csv output
# #######################################################################
 

print("# Writing to file")
                 
if settings.latex:
    
    def is_sig (v):
        if v <= 0.05:
            return "{\\bf %.4f}" % v
        else:
            return "%.4f" % v
            
    def is_zero (v):
        if v <= 1e-8:
            return "-"
        else:
            return "%.3f" % v
    
    print ("\\begin{tabular}{llllllllllll}\n\\hline")
    print ("Alignment & N & S & L & $AIC_c$ MH+S & \\multicolumn{3}{c}{$\\Delta AIC_c$} \\multicolumn{4}{c}{p-value}\\\\")
    print (" &  &  &  & & BUSTED & +S & +MH & +S+MH & BUSTED & +S & +MH\\\\\hline")
    
    sigByMethod  = Counter ()
    rankByMethod = Counter ()

    for file, result in sorted (by_file.items(), key=lambda record: record[1]['MH']['N'] * record[1]['MH']['S'] if 'MH' in record[1] else 0):
        if len (result) == 4:
            row = [file, str (result['MH']['N']), str (result['MH']['S']), "%.2f" % result['MH']['T']]
            if min ([k - result['MH']["AIC"] for k in [result['BUSTED']["AIC"],result['S']["AIC"], result['noS']["AIC"]]]) >= 0:
                row.append ("{\\bf %.1f}" % result['MH']["AIC"])
            else:
                row.append ("%.1f" % result['MH']["AIC"])
                
            cAIC = [['MH', result['MH']["AIC"]]]   
            
            for model in ['BUSTED','S','noS']:                
                row.append ("%.1f" % (result[model]["AIC"] - result['MH']["AIC"]))
                cAIC.append ([model,result[model]["AIC"]])
            
            best_model = sorted (cAIC, key=lambda record: record[1])[0]
            rankByMethod[best_model[0]] += 1
                
            
            for model in ['MH','BUSTED','S','noS']:
                row.append (is_sig(result[model]["p"]))
                if result[model]["p"] <= 0.05:
                    sigByMethod[model] += 1
                    
             
            print ('&'.join (row), "\\\\")
            
    print ("\\hline\n & & & & %d & %d & %d & %d & %d & %d & %d & %d  \\\\" % tuple ([rankByMethod [k] for k in ['MH','BUSTED','S','noS']] + [sigByMethod [k] for k in ['MH','BUSTED','S','noS']]))
    print ("\\hline\\end{tabular}")

    print ("\\begin{tabular}{llllllllllll}\n\\hline")
    print ("Alignment & $\\multicolumn{3}{c}{\\omega_3 (p_3)}$ & $CV (\\alpha)$ & $\\delta$ & $\\psi$ & \\multicolumn{3}{c}{Substitutuons (\\%)} & \\multicolumn{3}{c}{L} \\\\")
    print (" & +MH & +S &  & & & 1H & 2H & 3H & 1H & 2H & 3H\\\\\hline")

    for file, result in sorted (by_file.items(), key=lambda record: record[1]['MH']['N'] * record[1]['MH']['S'] if 'MH' in record[1] else 0):
        if len (result) == 4:
            row = [file]
            row.append ("%.3f (%.2f\\%%)" % (result['MH']["omega"]["omega"],100*result['MH']["omega"]["proportion"]))
            row.append ("%.3f (%.2f\\%%)" % (result['S']["omega"]["omega"],100*result['S']["omega"]["proportion"]))
            row.append ("%.3f" % result['MH']["SRV"])
            row.append (is_zero(result['MH']["delta"]))
            row.append (is_zero(result['MH']["psi"]))
            row.append ("%d (%.2f) & %d (%.2f) & %d (%.2f)" % (result['MH']["subs"][1], 100*result['MH']["subs"][1]/result['MH']["subs"][0],result['MH']["subs"][2], 100*result['MH']["subs"][2]/result['MH']["subs"][0],result['MH']["subs"][3], 100*result['MH']["subs"][3]/result['MH']["subs"][0]))
            row.append ("%.2f & %.2f & %.2f" % (result['MH']["L"][0],result['MH']["L"][1],result['MH']["L"][2]))
            print ('&'.join (row), "\\\\")
    print ("\\hline\\end{tabular}")

else: 
  
    headers = ['File','N','S','T','AIC-c','plain','+S','+MH','p','p_plain','p_s','p_MH']
    for model in ['BUSTED','S','noS','MH']:
        headers.append ("omega1_%s" % model)
        headers.append ("p1_%s" % model)
        headers.append ("omega2_%s" % model)
        headers.append ("p2_%s" % model)
        headers.append ("omega3_%s" % model)
        headers.append ("p3_%s" % model)
    #end for

    for model in ['S','MH']:
         headers.append ("SRV_%s" % model)
    #end for
   
    headers.extend (['delta','psi', 'S1', 'S2', 'S3', 'SP1' , 'SP2', 'SP3', 'L1', 'L2','L3','EB'])

    # Clear output, print headers
    with open(settings.output, 'w', newline='') as csvfile:
        output_writer = csv.writer (csvfile)
        output_writer.writerow (headers)
    #end for

    pv = Counter ()

    for file, result in sorted (by_file.items(), key=lambda record: record[1]['MH']['N'] * record[1]['MH']['S'] if 'MH' in record[1] else 0):
        if len (result) == 4:
            row = [file, str (result['MH']['N']), str (result['MH']['S']), "%.2f" % result['MH']['T']]
            row.append ("%.5f" % result['MH']["AIC"])
            row.append ("%.5f" % (result['BUSTED']["AIC"] - result['MH']["AIC"]))
            row.append ("%.5f" % (result['S']["AIC"] - result['MH']["AIC"]))
            row.append ("%.5f" % (result['noS']["AIC"] - result['MH']["AIC"]))
            row.append ("%.5f" % result['MH']["p"])      # p = BUSTEDS-MH
            row.append ("%.5f" % result['BUSTED']["p"])  # p_plain = BUSTED
            row.append ("%.5f" % result['S']["p"])       # p_s = BUSTEDS
            row.append ("%.5f" % result['noS']["p"])     # p_MH = BUSTED-MH
            for model in ['BUSTED','S','noS','MH']:
                row.append ("%.5f" % (result[model]["omega1"]["omega"]))
                row.append ("%.5f" % (result[model]["omega1"]["proportion"]))
                row.append ("%.5f" % (result[model]["omega2"]["omega"]))
                row.append ("%.5f" % (result[model]["omega2"]["proportion"]))
                row.append ("%.5f" % (result[model]["omega"]["omega"]))
                row.append ("%.5f" % (result[model]["omega"]["proportion"]))

 
            for model in ['S','MH']:
                row.append ("%.5f" % (result[model]["SRV"]))
            #end for
  
            row.append ("%.5f" % result['MH']["delta"])
            row.append ("%.5f" % result['MH']["psi"])
            row.append ("%d" % (result['MH']["subs"][1]))
            row.append ("%d" % (result['MH']["subs"][2]))
            row.append ("%d" % (result['MH']["subs"][3]))
            row.append ("%g" % (result['MH']["subs"][1]/result['MH']["subs"][0]))
            row.append ("%g" % (result['MH']["subs"][2]/result['MH']["subs"][0]))
            row.append ("%g" % (result['MH']["subs"][3]/result['MH']["subs"][0]))
            row.append ("%.5f" % (result['MH']["L"][0]))
            row.append ("%.5f" % (result['MH']["L"][1]))
            row.append ("%.5f" % (result['MH']["L"][2]))
            row.append ("%.d" % result['MH']["EB"])
            for model in ['BUSTED','S','noS','MH']:
                pv[model] += 1 if result [model]['p'] <= 0.05 else 0
            pv['N'] += 1
            #print (result['MH']['p'], file = sys.stderr)

            with open(settings.output, 'a+', newline='') as csvfile:
                output_writer = csv.writer (csvfile)
                #output_writer.writerow (headers)
                #output_writer = csv.writer (sys.stdout)
                output_writer.writerow (row)
            #end with
        
    print (pv, file = sys.stderr)

# #######################################################################
# End of file
# #######################################################################
                     
