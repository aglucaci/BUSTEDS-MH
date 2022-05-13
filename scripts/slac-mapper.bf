// Load libraries
LoadFunctionLibrary ("libv3/tasks/trees.bf");
LoadFunctionLibrary ("libv3/IOFunctions.bf");

// Read the SLAC json
// First command line argument
slac_json = io.ReadFromOrCreate ("Load SLAC fit file", {})["value"];

site_count = (slac_json["input"])["number of sites"];
tree_str = ((slac_json["input"])["trees"])[0];

// Load topology from input
Topology T = tree_str;

// store node order
branch_names = BranchName (T,-1);
branch_names_map = {};

// Get branch attributes
slac_json = (slac_json["branch attributes"])[0];
site_counts = {};
bl = {};

for (i; in; branch_names) {
    branch_names_map + i;
    bl[i] = (slac_json[i])["Global MG94xREV"]; 
}

parent_index      = Abs (T);
subs_by_count     = {};
sites_by_count    = {};
branches_by_count = {}; 

// Main logic here
for (s = 0; s < site_count; s+=1) {
    site_map = {};
    subs = {{0,0,0,0}};
    site_subs = {};
    for (i,b; in; branch_names)    {
        p = parent_index[i];
        is_i = 2*((b[0][3] && 1)=="NODE");
        my_state = ((slac_json[b])["codon"])[0][s];
        NSS = ((slac_json[b])["synonymous substitution count"])[0][s] + ((slac_json[b])["nonsynonymous substitution count"])[0][s];
        if (NSS > 0) {
            subs_by_count [NSS] += 1;
            if (sites_by_count/NSS == FALSE) {
                sites_by_count [NSS] = {};
            }
        
            (sites_by_count [NSS])[s] = 1;
            if (branches_by_count/NSS == FALSE) {
                branches_by_count [NSS] = {};
            }
            (branches_by_count[NSS])[b] = 1;
        }
    } 
    
}

// Create empty output file
json_out = {};

for (subs, count; in; subs_by_count) {
    json_out [subs] = {
        'count' : count,
        'sites' : Abs (sites_by_count [subs]),
        'branches' : Abs (branches_by_count [subs])
    };
}


// Output file is second command line argument
io.SpoolJSON (json_out, io.PromptUserForFilePath("Save JSON to"));


// End of file