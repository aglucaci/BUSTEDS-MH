"""
Combine analysis results

Author:
    Sergei L Kosakovsky Pond (spond@temple.edu)

Version:
    v0.0.1 (2021-01-17)


"""

import argparse
import csv
import random
import os
import json
import sys
import re
import math
import numpy
from scipy.stats import chi2

import progressbar

bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)

from collections import Counter

random.seed ()

arguments = argparse.ArgumentParser(description='Combine alignments into a single file, adding a reference sequence as well')

arguments.add_argument('-i', '--input',  help = 'Directories ', required = True, type = str, nargs = '*')
arguments.add_argument('-s', '--simulations',  help = 'Directories ', action = 'store_true' )
arguments.add_argument('-l', '--latex',  help = 'Make LATeX output ', action = 'store_true' )

settings = arguments.parse_args()

by_file = {}

timer = 0
count = 0
tags = {}

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
    
nucs = set (['A','C','G','T'])
    
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
 
 
i = 0

    
for basedir in settings.input:
     for root, dirs, files in os.walk(basedir):
        for each_file in files:
            file_name, file_ext =  os.path.splitext (each_file)
            #print (file_ext, file = sys.stderr)
            if file_ext == '.json':
                parts = file_name.split ('.')
                #print (parts)
                with open (os.path.join (root, each_file), "r") as fh:
                    try:
                        i += 1
                        bar.update(i)

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
 
                    
    
         
    
omegas = {'MH' : [], 'S' : []}
weights = {'MH' : [], 'S' : []}
deltas = []
psi = []
pvalues = {'MH' : [], 'S' : [], 'MA' : []}
aic = {'MH' : [], 'S' : []}
lrt_pv = []


def p_rate (pv, cutoff = 0.05):
    if len (pv) == 0: return "N/A"
    return "%.2g" % (len ([p for p in pv if p <= cutoff]) / float (len (pv)))

def iqr (rates, cutoff = 0.05):
    iqr_v = numpy.quantile (rates, [0.25,0.75])
    return "%.2f-%.2f" % (iqr_v[0], iqr_v[1])

for file, result in sorted (by_file.items(), key=lambda record: record[1]['MH']['N'] * record[1]['MH']['S'] if 'MH' in record[1] else 0):
    if len (result) >= 2:
        
        min_AIC = 1e100
        
        for model in ['MH','S']:                
            omegas[model].append (result[model]['omega']['omega'])
            weights[model].append (result[model]['omega']['proportion'])
            pvalues[model].append (result[model]['p'])
            aic[model].append (result[model]['AIC'])
            min_AIC = min (min_AIC,result[model]['AIC'])
            
        rw_aic = []   
        for model in ['MH','S']:
            rw_aic.append (math.exp ((min_AIC-result[model]['AIC']) * 0.5))
            
        lrt_pv.append (1-chi2.cdf (max (0,2*(result['MH']['logL'] - result['S']['logL'])),2))
            
        rw_sum = sum (rw_aic)
        rw_aic = [k / rw_sum for k in rw_aic]
        rw_p = result['MH']['p'] * rw_aic [0] + result['S']['p'] * rw_aic [1]
        pvalues['MA'].append (rw_p)
        if rw_p > 0.05:
            print (rw_aic, result['MH']['p'], result['S']['p'],rw_p)
                
        deltas.append (result['MH']['delta'])
        psi.append (result['MH']['psi'])
        
mh_better = [i for i,p in enumerate (aic["MH"]) if p < aic["S"][i]]
s_better = [i for i,p in enumerate (aic["MH"]) if p >= aic["S"][i]]
p_mh = [pvalues['MH'][i] for i in mh_better]
p_s = [pvalues['S'][i] for i in s_better]



print ("omega (weight) & $ %s $  ($%.2g\\%%$)& $ %s $  ($%.2g\\%%$) & delta & $ %s $ & psi & $ %s $  & %s (%d) & %s (%d) & %s & %s"  % (iqr (omegas['S']), 100*numpy.mean (weights['S']), iqr (omegas['MH']), 100*numpy.mean (weights['MH']),iqr (deltas), iqr (psi), p_rate (pvalues['S']), len (p_s), p_rate (pvalues['MH']), len (p_mh), p_rate (pvalues['MA']), p_rate (lrt_pv)))
       
print ("\n\n", len (psi))     
    

                        
