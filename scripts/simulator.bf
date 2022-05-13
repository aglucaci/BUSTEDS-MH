LoadFunctionLibrary("libv3/all-terms.bf");
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/models/parameters.bf");

KeywordArgument  ("fit", "Load the BUSTED-MH fit file");
sim.filepath = io.LoadFile("Load the BUSTED-MH fit file");

// print to user
console.log (">Loaded fit file from `sim.filepath`");

// list of likelihood function IDs
sim.lf = (utility.GetListOfLoadedLikelihoodFunctions (None));

// assert statement, checks the array.
// how is this stored?
assert (utility.Array1D(sim.lf) == 1, "The fit file must contain a single likelihood function");

// the ID of the single loaded likelihood function
sim.lf  = sim.lf [0];  

// get the list of likelihood function parameters
GetString (sim.parameters, ^sim.lf,-1);

/* example
            {
         "Base frequencies":  {
        {"busted.test_pi"} 
          },
         "Categories":  {
        {"busted.test.rv_gdd"} 
          },
         "Compute Template":"",
         "Datafilters":  {
        {"busted.filter.default"} 
          },
         "Global Constrained":  {
        {"busted.test.rv_gdd_norm", "busted.test.theta_AG"} 
          },
         "Global Independent":  {
        {"busted.test.theta_AC", "busted.test.omega1", "busted.test.theta_AT", "busted.test.theta_CG", "busted.test.psi_islands", "busted.test.omega3", "busted.test.delta", "busted.test.psi", "busted.test.omega2", "busted.test.theta_GT", "busted.test.rv_gdd_weights_1", "busted.test.bsrel_mixture_aux_1", "busted.test.rv_gdd_rates_0", "busted.test.theta_CT", "busted.test.rv_gdd_weights_0", "busted.test.rv_gdd_rates_2", "busted.test.rv_gdd_rates_1", "busted.test.bsrel_mixture_aux_0"} 
          },
         "Local Constrained":  {
          },
         "Local Independent":  {
        {"iOBFdqGX.tree_id_0.YAK.t", "iOBFdqGX.tree_id_0.MEL.t", "iOBFdqGX.tree_id_0.SIM.t", "iOBFdqGX.tree_id_0.MA.t", "iOBFdqGX.tree_id_0.Node7.t", "iOBFdqGX.tree_id_0.Node5.t", "iOBFdqGX.tree_id_0.Node3.t", "iOBFdqGX.tree_id_0.ERE.t", "iOBFdqGX.tree_id_0.Node2.t", "iOBFdqGX.tree_id_0.PSE.t", "iOBFdqGX.tree_id_0.PS.t", "iOBFdqGX.tree_id_0.Node14.t", "iOBFdqGX.tree_id_0.PER.t", "iOBFdqGX.tree_id_0.Node13.t", "iOBFdqGX.tree_id_0.MIR.t", "iOBFdqGX.tree_id_0.Node12.t", "iOBFdqGX.tree_id_0.SUB.t", "iOBFdqGX.tree_id_0.AMB.t", "iOBFdqGX.tree_id_0.Node19.t", "iOBFdqGX.tree_id_0.Node11.t", "iOBFdqGX.tree_id_0.Node1.t", "iOBFdqGX.tree_id_0.MET.t", "iOBFdqGX.tree_id_0.CRA.t", "iOBFdqGX.tree_id_0.NIG.t", "iOBFdqGX.tree_id_0.MIM.t", "iOBFdqGX.tree_id_0.ADI.t", "iOBFdqGX.tree_id_0.PIC.t", "iOBFdqGX.tree_id_0.SIL.t", "iOBFdqGX.tree_id_0.HET.t", "iOBFdqGX.tree_id_0.Node36.t", "iOBFdqGX.tree_id_0.PLA.t", "iOBFdqGX.tree_id_0.DIF.t", "iOBFdqGX.tree_id_0.Node39.t", "iOBFdqGX.tree_id_0.Node35.t", "iOBFdqGX.tree_id_0.AFF.t", "iOBFdqGX.tree_id_0.Node34.t", "iOBFdqGX.tree_id_0.Node32.t", "iOBFdqGX.tree_id_0.Node30.t", "iOBFdqGX.tree_id_0.Node28.t", "iOBFdqGX.tree_id_0.Node26.t", "iOBFdqGX.tree_id_0.Node24.t", "iOBFdqGX.tree_id_0.Node22.t", "iOBFdqGX.tree_id_0.LEB.t"} 
          },
         "Models":  {
        {"busted.test"} 
          },
         "Trees":  {
        {"iOBFdqGX.tree_id_0"} 
          }
        }
*/

// Get the global parameters, w, srv, dh, th, th_si, nt rates, substition params, etc.
sim.globals = sim.parameters [terms.parameters.global_independent];
console.log (">The likelihood function `sim.lf` has `utility.Array1D(sim.globals)` global parameters");
//console.log ("\nGlobal parameters and their MLEs (ranges)");

// Enter for loop, looping over global parameters.
for (sim.param; in; sim.globals) {
    sim.param.MLE = Eval (sim.param);
    sim.param.range = parameters.GetRange(sim.param);
    //console.log ("\t" + sim.param + " = " + sim.param.MLE + " [" + sim.param.range [terms.lower_bound] + ", " +sim.param.range [terms.upper_bound] + "]");
    
    ExecuteCommands ('
        KeywordArgument  ("`sim.param`", "Value for `sim_param` to use during simulations", `sim.param.MLE`);
        sim.use_this_value = io.PromptUser ("Value for `sim.param` to use during simulations",
                                            `sim.param.MLE`,
                                            `sim.param.range [terms.lower_bound]`,
                                            `sim.param.range [terms.upper_bound]`,
                                            0);
    ');
    
    ^sim.param = sim.use_this_value;
}

// How many replicates to make.
KeywordArgument  ("replicates", "How many replicates to generate", 10);
sim.replicates = io.PromptUser ("How many replicates to generate", 10, 1, 100000, TRUE);
console.log ("Will generate `sim.replicates` replicates");

// Output file prefix.
KeywordArgument ("output", "Write the function used for data simulation and replicates to this path");
sim.path = io.PromptUserForFilePath ( "Write the function used for data simulation and replicates to this path");

Export (sim.lf_export, ^sim.lf);

fprintf (sim.path, CLEAR_FILE, sim.lf_export);

for (i = 1; i <= sim.replicates; i+=1) {
    io.ReportProgressBar("", "Generating replicate " + i + "/" +  sim.replicates);
    // Heavy lifting here.
    DataSet sim.simulated = SimulateDataSet (^sim.lf);

    // What does the filter do?
    DataSetFilter sim.simulated.filter = CreateFilter (sim.simulated,1);
  
    // Report to user.
    fprintf (sim.path + ".replicate." + i, CLEAR_FILE, sim.simulated.filter);
}

io.ClearProgressBar();


// End of file.

