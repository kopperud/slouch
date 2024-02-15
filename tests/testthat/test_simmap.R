context("testing simmap likelihood")

## hard-code a SIMMAP tree
phy <- structure(list(edge = structure(c(44L, 44L, 45L, 46L, 47L, 48L, 
                                         49L, 50L, 51L, 52L, 53L, 54L, 54L, 55L, 55L, 53L, 56L, 56L, 52L, 
                                         57L, 57L, 58L, 58L, 51L, 50L, 59L, 60L, 61L, 62L, 62L, 61L, 60L, 
                                         59L, 49L, 63L, 64L, 65L, 66L, 67L, 68L, 68L, 67L, 69L, 69L, 70L, 
                                         71L, 71L, 72L, 72L, 70L, 66L, 65L, 64L, 73L, 73L, 63L, 74L, 75L, 
                                         75L, 74L, 76L, 76L, 48L, 77L, 77L, 47L, 78L, 78L, 79L, 80L, 81L, 
                                         82L, 82L, 81L, 83L, 84L, 84L, 83L, 80L, 79L, 46L, 45L, 85L, 85L, 
                                         1L, 45L, 46L, 47L, 48L, 49L, 50L, 51L, 52L, 53L, 54L, 2L, 55L, 
                                         23L, 24L, 56L, 14L, 15L, 57L, 4L, 58L, 10L, 11L, 30L, 59L, 60L, 
                                         61L, 62L, 7L, 8L, 32L, 26L, 22L, 63L, 64L, 65L, 66L, 67L, 68L, 
                                         5L, 18L, 69L, 6L, 70L, 71L, 12L, 72L, 19L, 20L, 13L, 31L, 25L, 
                                         73L, 9L, 27L, 74L, 75L, 16L, 17L, 76L, 28L, 29L, 77L, 3L, 21L, 
                                         78L, 33L, 79L, 80L, 81L, 82L, 34L, 40L, 83L, 84L, 36L, 39L, 38L, 
                                         35L, 37L, 41L, 85L, 42L, 43L), dim = c(84L, 2L)), edge.length = c(27.1899206830978, 
                                                                                                           3.27516624544859, 1.24340673072568, 3.14305207598277, 1.6492431943059, 
                                                                                                           0.918495305532218, 0.676900710988054, 2.31360760892056, 0.781373936992214, 
                                                                                                           6.29905043254729, 3.53564415827989, 3.3539802688241, 1.50538719616652, 
                                                                                                           1.84859307269838, 1.84859307293679, 1.30896923911869, 5.58065518739224, 
                                                                                                           5.58065518739224, 6.78370486610003, 6.40496999832988, 4.98637419116497, 
                                                                                                           1.41859580630064, 1.41859580630064, 13.9700488012952, 1.15317506787492, 
                                                                                                           4.9422154805541, 2.69835678109527, 5.88634803197794, 1.60356104657054, 
                                                                                                           1.60356104724407, 7.48990907708257, 10.188265863359, 15.1304813441277, 
                                                                                                           0.296001798122351, 2.66731694515798, 0.874660799215309, 2.35726635359453, 
                                                                                                           0.680342813456981, 1.8211930157384, 8.26377539496422, 8.26377539496422, 
                                                                                                           2.87572802776304, 7.2092403859362, 0.698727641859382, 1.07334864132804, 
                                                                                                           5.43716409749985, 3.19540856140852, 2.2417555366978, 2.24175553792715, 
                                                                                                           6.51051274144564, 10.7653112264007, 13.1225775730603, 1.36101562911711, 
                                                                                                           12.6362227455974, 12.6362227487564, 9.67617326985942, 4.25246768498421, 
                                                                                                           2.73591436501741, 2.73591436501741, 2.98354652861357, 4.00483552007675, 
                                                                                                           4.00483552007675, 2.01095457867384, 15.868097840786, 15.8680978440762, 
                                                                                                           4.39560088924909, 15.1326947284911, 7.80886816276765, 0.899028033828735, 
                                                                                                           1.53370201078057, 0.527478770575071, 4.36361775184371, 4.36361775116236, 
                                                                                                           1.00235363149419, 0.536495257891715, 3.35224763460159, 3.35224763460159, 
                                                                                                           3.88874289157391, 6.42479853405952, 7.32382656793594, 22.6713477031208, 
                                                                                                           7.93035265209675, 15.9844017724037, 15.9844017724037), Nnode = 42L, 
                      tip.label = c("Antilocapra_americana", "Addax_nasomaculatus", 
                                    "Aepyceros_melampus", "Alcelaphus_buselaphus_buselaphus", 
                                    "Antidorcas_marsupialis", "Antilope_cervicapra", "Cephalophus_natalensis", 
                                    "Cephalophus_nigrifrons", "Madoqua_kirkii", "Connochaetes_gnou", 
                                    "Connochaetes_taurinus_1", "Eudorcas_rufifrons_2", "Gazella_dorcas_pelzelnii", 
                                    "Hippotragus_equinus", "Hippotragus_niger", "Kobus_ellipsiprymnus", 
                                    "Kobus_leche", "Litocranius_walleri", "Nanger_dama", "Nanger_granti", 
                                    "Nesotragus_moschatus", "Oreotragus_oreotragus", "Oryx_dammah", 
                                    "Oryx_gazella", "Ourebia_ourebi", "Philantomba_monticola_1", 
                                    "Rahicerus_campestris", "Redunca_arundinum", "Redunca_fulvorufula", 
                                    "Rupicapra_rupicapra_FJ207538", "Saiga_tatarica", "Sylvicapra_grimmia", 
                                    "Syncerus_caffer", "Tragelaphus_oryx", "Tragelaphus_angasii", 
                                    "Tragelaphus_eurycerus", "Tragelaphus_imberbis", "Tragelaphus_scriptus_3", 
                                    "Tragelaphus_spekii", "Tragelaphus_strepsiceros", "Capreolus_capreolus", 
                                    "Giraffa_camelopardalis_angolensis_NC012100", "Okapia_johnstoni"
                      ), node.label = c("MF", "MF", "MF", "MF", "MF", "MF", "MF", 
                                        "MF", "MF", "MF", "MF", "MF", "Gr", "Gr", "Gr", "MF", "Br", 
                                        "Br", "Br", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", 
                                        "MF", "MF", "MF", "MF", "MF", "Gr", "MF", "MF", "MF", "MF", 
                                        "MF", "MF", "MF", "MF", "Br"), maps = list(c(MF = 22.9800362350042, 
                                                                                     Gr = 1.75707098691584, MF = 2.45281346117775), c(MF = 0.5040779628898, 
                                                                                                                                      Br = 2.77108828255879), c(Br = 1.24340673072568), c(Br = 0.796910970538387, 
                                                                                                                                                                                          MF = 2.34614110544438), c(MF = 1.6492431943059), c(MF = 0.918495305532218), 
                                                                                   c(MF = 0.676900710988054), c(MF = 2.31360760892056), 
                                                                                   c(MF = 0.727728042675122, Gr = 0.0536458943170928), c(Gr = 6.29905043254729), 
                                                                                   c(Gr = 2.35212643509007, MF = 1.18351772318982), c(MF = 3.3539802688241), 
                                                                                   c(MF = 1.50538719616652), c(MF = 1.84859307269838), c(MF = 1.84859307293679), 
                                                                                   c(Gr = 1.30896923911869), c(Gr = 5.58065518739224), c(Gr = 5.58065518739224), 
                                                                                   c(Gr = 6.78370486610003), c(Gr = 6.40496999832988), c(Gr = 4.98637419116497), 
                                                                                   c(Gr = 1.41859580630064), c(Gr = 1.41859580630064), c(MF = 13.9700488012952), 
                                                                                   c(MF = 1.15317506787492), c(MF = 4.08117839167061, Br = 0.86103708888349
                                                                                   ), c(Br = 2.69835678109527), c(Br = 5.88634803197794), 
                                                                                   c(Br = 1.60356104657054), c(Br = 1.60356104724407), c(Br = 7.48990907708257), 
                                                                                   c(Br = 10.188265863359), c(MF = 5.5063469667341, Br = 3.65581881765503, 
                                                                                                              MF = 5.96831555973852), c(MF = 0.296001798122351), c(MF = 2.66731694515798), 
                                                                                   c(MF = 0.874660799215309), c(MF = 2.35726635359453), 
                                                                                   c(MF = 0.680342813456981), c(MF = 1.8211930157384), c(MF = 8.26377539496422), 
                                                                                   c(MF = 1.61929623891263, Br = 6.64447915605159), c(MF = 2.87572802776304), 
                                                                                   c(MF = 7.2092403859362), c(MF = 0.698727641859382), c(MF = 1.07334864132804), 
                                                                                   c(MF = 5.43716409749985), c(MF = 3.19540856140852), c(MF = 2.2417555366978), 
                                                                                   c(MF = 2.24175553792715), c(MF = 6.51051274144564), c(MF = 10.7653112264007), 
                                                                                   c(MF = 13.1225775730603), c(MF = 1.36101562911711), c(MF = 12.0833300497328, 
                                                                                                                                         Br = 0.552892695864521), c(MF = 12.6362227487564), c(MF = 9.67617326985942), 
                                                                                   c(MF = 4.25246768498421), c(MF = 2.73591436501741), c(MF = 2.73591436501741), 
                                                                                   c(MF = 0.417969927947861, Gr = 2.56557660066571), c(Gr = 4.00483552007675), 
                                                                                   c(Gr = 4.00483552007675), c(MF = 1.29650540940574, Br = 0.714449169268102
                                                                                   ), c(Br = 12.9977917003104, MF = 2.87030614047554), c(Br = 15.8680978440762), 
                                                                                   c(MF = 4.39560088924909), c(MF = 6.49233349799687, Gr = 8.64036123049427
                                                                                   ), c(MF = 0.844801102347294, Br = 0.861372853802622, 
                                                                                        MF = 6.10269420661773), c(MF = 0.899028033828735), c(MF = 1.53370201078057), 
                                                                                   c(MF = 0.527478770575071), c(MF = 4.36361775184371), 
                                                                                   c(MF = 1.22905756392906, Br = 3.1345601872333), c(MF = 1.00235363149419), 
                                                                                   c(MF = 0.536495257891715), c(MF = 3.35224763460159), 
                                                                                   c(MF = 3.35224763460159), c(MF = 3.00711901712153, Br = 0.881623874452379
                                                                                   ), c(MF = 0.34037749717795, Br = 6.08442103688157), c(MF = 7.32382656793594), 
                                                                                   c(Br = 22.6713477031208), c(Br = 7.93035265209675), c(Br = 15.9844017724037), 
                                                                                   c(Br = 0.55592172339438, MF = 6.44902649182669, Br = 8.97945355718265
                                                                                   )), mapped.edge = structure(c(0, 2.77108828255879, 1.24340673072568, 
                                                                                                                 0.796910970538387, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                 0, 0, 0, 0, 0, 0, 0, 0, 0.86103708888349, 2.69835678109527, 
                                                                                                                 5.88634803197794, 1.60356104657054, 1.60356104724407, 7.48990907708257, 
                                                                                                                 10.188265863359, 3.65581881765503, 0, 0, 0, 0, 0, 0, 0, 6.64447915605159, 
                                                                                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.552892695864521, 0, 
                                                                                                                 0, 0, 0, 0, 0, 0, 0, 0.714449169268102, 12.9977917003104, 
                                                                                                                 15.8680978440762, 0, 0, 0.861372853802622, 0, 0, 0, 0, 3.1345601872333, 
                                                                                                                 0, 0, 0, 0, 0.881623874452379, 6.08442103688157, 0, 22.6713477031208, 
                                                                                                                 7.93035265209675, 15.9844017724037, 9.53537528057703, 1.75707098691584, 
                                                                                                                 0, 0, 0, 0, 0, 0, 0, 0.0536458943170928, 6.29905043254729, 
                                                                                                                 2.35212643509007, 0, 0, 0, 0, 1.30896923911869, 5.58065518739224, 
                                                                                                                 5.58065518739224, 6.78370486610003, 6.40496999832988, 4.98637419116497, 
                                                                                                                 1.41859580630064, 1.41859580630064, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.56557660066571, 4.00483552007675, 
                                                                                                                 4.00483552007675, 0, 0, 0, 0, 8.64036123049427, 0, 0, 0, 
                                                                                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25.432849696182, 
                                                                                                                 0.5040779628898, 0, 2.34614110544438, 1.6492431943059, 0.918495305532218, 
                                                                                                                 0.676900710988054, 2.31360760892056, 0.727728042675122, 0, 
                                                                                                                 1.18351772318982, 3.3539802688241, 1.50538719616652, 1.84859307269838, 
                                                                                                                 1.84859307293679, 0, 0, 0, 0, 0, 0, 0, 0, 13.9700488012952, 
                                                                                                                 1.15317506787492, 4.08117839167061, 0, 0, 0, 0, 0, 0, 11.4746625264726, 
                                                                                                                 0.296001798122351, 2.66731694515798, 0.874660799215309, 2.35726635359453, 
                                                                                                                 0.680342813456981, 1.8211930157384, 8.26377539496422, 1.61929623891263, 
                                                                                                                 2.87572802776304, 7.2092403859362, 0.698727641859382, 1.07334864132804, 
                                                                                                                 5.43716409749985, 3.19540856140852, 2.2417555366978, 2.24175553792715, 
                                                                                                                 6.51051274144564, 10.7653112264007, 13.1225775730603, 1.36101562911711, 
                                                                                                                 12.0833300497328, 12.6362227487564, 9.67617326985942, 4.25246768498421, 
                                                                                                                 2.73591436501741, 2.73591436501741, 0.417969927947861, 0, 
                                                                                                                 0, 1.29650540940574, 2.87030614047554, 0, 4.39560088924909, 
                                                                                                                 6.49233349799687, 6.94749530896502, 0.899028033828735, 1.53370201078057, 
                                                                                                                 0.527478770575071, 4.36361775184371, 1.22905756392906, 1.00235363149419, 
                                                                                                                 0.536495257891715, 3.35224763460159, 3.35224763460159, 3.00711901712153, 
                                                                                                                 0.34037749717795, 7.32382656793594, 0, 0, 0, 6.44902649182669
                                                                                   ), dim = c(84L, 3L), dimnames = list(c("44,1", "44,45", "45,46", 
                                                                                                                          "46,47", "47,48", "48,49", "49,50", "50,51", "51,52", "52,53", 
                                                                                                                          "53,54", "54,2", "54,55", "55,23", "55,24", "53,56", "56,14", 
                                                                                                                          "56,15", "52,57", "57,4", "57,58", "58,10", "58,11", "51,30", 
                                                                                                                          "50,59", "59,60", "60,61", "61,62", "62,7", "62,8", "61,32", 
                                                                                                                          "60,26", "59,22", "49,63", "63,64", "64,65", "65,66", "66,67", 
                                                                                                                          "67,68", "68,5", "68,18", "67,69", "69,6", "69,70", "70,71", 
                                                                                                                          "71,12", "71,72", "72,19", "72,20", "70,13", "66,31", "65,25", 
                                                                                                                          "64,73", "73,9", "73,27", "63,74", "74,75", "75,16", "75,17", 
                                                                                                                          "74,76", "76,28", "76,29", "48,77", "77,3", "77,21", "47,78", 
                                                                                                                          "78,33", "78,79", "79,80", "80,81", "81,82", "82,34", "82,40", 
                                                                                                                          "81,83", "83,84", "84,36", "84,39", "83,38", "80,35", "79,37", 
                                                                                                                          "46,41", "45,85", "85,42", "85,43"), c("Br", "Gr", "MF"))), 
                      Q = structure(c(-0.0398636047637968, 0, 0.0398636047637968, 
                                      0, -0.0179715898466568, 0.0179715898466568, 0.0398636047637968, 
                                      0.0179715898466568, -0.0578351946104536), dim = c(3L, 3L), dimnames = list(
                                        c("Br", "Gr", "MF"), c("Br", "Gr", "MF"))), logL = structure(-36.797331003745, df = 3L)), class = c("simmap", 
                                                                                                                                            "phylo"), order = "cladewise", map.order = "right-to-left")

### alternative likelihood implementation

sd_postorder <- function(node_index, edge, tree, continuousChar,
                         mu, V, log_norm_factor, subedges_lengths, sigma_squared, alpha, theta){
  ntip = length(tree$tip.label)
  
  # if is internal node
  if (node_index > ntip){
    
    left_edge  = which(edge[,1] == node_index)[1] # index of left child edge
    right_edge = which(edge[,1] == node_index)[2] # index of right child edge
    left = edge[left_edge,2] # index of left child node
    right = edge[right_edge,2] # index of right child node
    
    output_left <- sd_postorder(left, edge, tree, continuousChar,
                                mu, V, log_norm_factor, subedges_lengths, sigma_squared, alpha, theta)
    mu <- output_left[[1]]
    V <- output_left[[2]]
    log_norm_factor <- output_left[[3]]
    
    output_right <- sd_postorder(right, edge, tree, continuousChar,
                                 mu, V, log_norm_factor, subedges_lengths, sigma_squared, alpha, theta)
    mu <- output_right[[1]]
    V <- output_right[[2]]
    log_norm_factor <- output_right[[3]]
    
    sub_bl_left = subedges_lengths[left_edge][[1]] # all subedges of left child edge
    sub_bl_right = subedges_lengths[right_edge][[1]] # all subedges of right child edge
    
    # for the sake of readability, computation of variance, mean, and log_nf are done in separate loops
    # 1) variance of the normal variable: this branch (v_left) and the subtree (V[left])
    ## Is 'delta_left* exp(2.0 * alpha * bl_left)' added in each sub-edge?
    
    delta_left = V[left]
    v_left = 0 # initialise v_left
    for (i in rev(1:length(sub_bl_left))){
      state <- names(sub_bl_left)[i]
      delta_t <- sub_bl_left[[i]]
      v_left = sigma_squared/(2*alpha) * expm1(2.0*alpha * delta_t) 
      delta_left = v_left + delta_left * exp(2.0 * alpha * delta_t)
    }
    
    delta_right = V[right]
    v_right = 0 # initialise v_right
    for (i in rev(1:length(sub_bl_right))){
      state <- names(sub_bl_right)[i]
      v_right = sigma_squared/(2*alpha) *expm1(2.0*alpha*sub_bl_right[[i]])
      delta_right = v_right + delta_right * exp(2.0 * alpha * sub_bl_right[[i]])
    }
    
    var_left = delta_left
    var_right = delta_right
    
    # 2) mean of the normal variable
    mean_left = mu[left]
    for (i in rev(1:length(sub_bl_left))){
      state <- names(sub_bl_left)[i]
      mean_left = exp(alpha*sub_bl_left[[i]])*(mean_left - theta[[state]]) + theta[[state]]
    }
    
    mean_right = mu[right]
    for (i in rev(1:length(sub_bl_right))){
      state <- names(sub_bl_right)[i]
      mean_right = exp(alpha*sub_bl_right[[i]])*(mean_right - theta[[state]]) + theta[[state]]
    }
    
    ## compute the mean and variance of the node
    mean_ancestor = (mean_left * var_right + mean_right * var_left) / (var_left + var_right)
    mu[node_index] = mean_ancestor
    var_node = (var_left * var_right) / (var_left + var_right)
    V[node_index] = var_node
    
    ## compute the normalizing factor, the left-hand side of the pdf of the normal variable
    log_nf_left = 0
    for (i in rev(1:length(sub_bl_left))){
      state <- names(sub_bl_left)[i]
      delta_t <- sub_bl_left[[i]]
      log_nf_left = log_nf_left + delta_t * alpha
    }
    
    log_nf_right = 0
    for (i in rev(1:length(sub_bl_right))){
      state <- names(sub_bl_right)[i]
      log_nf_right = log_nf_right + sub_bl_right[[i]] * alpha
    }
    
    contrast = mean_left - mean_right
    a = -(contrast*contrast / (2*(var_left+var_right)))
    b = log(2*pi*(var_left+var_right))/2.0
    log_nf = log_nf_left + log_nf_right + a - b
    log_norm_factor[node_index] = log_nf
    
    return(list(mu, V, log_norm_factor))
  }
  # if is tip
  else{
    species = tree$tip.label[node_index]
    
    mu[node_index] = continuousChar[[which(names(continuousChar) == species)]]
    V[node_index] = 0.0 ## if there is no observation error
    
    return(list(mu, V, log_norm_factor))
  }
}

sd_logL_pruning <- function(tree, continuousChar, sigma_squared, alpha, theta){
  ntip = length(tree$tip.label) # number of tips
  edge = tree$edge # equals tree[:edge] in Julia
  n_edges = length(edge[,1]) # number of edges
  max_node_index = max(tree$edge) # total number of nodes
  
  V = numeric(max_node_index)
  mu = numeric(max_node_index)
  log_norm_factor = numeric(max_node_index)
  
  subedges_lengths = tree$maps
  
  root_index = ntip + 1
  
  output <- sd_postorder(root_index, edge, tree, continuousChar,
                         mu, V, log_norm_factor, subedges_lengths, sigma_squared, alpha, theta)
  mu <- output[[1]]
  V <- output[[2]]
  log_norm_factor <- output[[3]]
  
  ## assume root value equal to theta
  mu_root = mu[root_index]
  v_root = V[root_index]
  left_edge_from_root <- which(edge[,1] == ntip+1)[1] # obtain left child edge index of root node
  left_subedges_from_root <- subedges_lengths[[left_edge_from_root]] # obtain sub-edge lengths
  root_state = names(tail(left_subedges_from_root))[[1]] # obtain root state, assuming it equals last state at left child edge
  lnl = dnorm(theta[[root_state]], mean = mu_root, sd = sqrt(v_root), log = TRUE)
  
  ## add norm factor
  for (log_nf in log_norm_factor){
    lnl = lnl + log_nf
  }
  return(lnl)
}

# test with slouch data set
data("artiodactyla")
data("neocortex")

# convert continuous data to read.nexus.data() format
brain <- list()
for (i in 1:length(neocortex$brain_mass_g_log_mean)){
  sp <- neocortex$species[i]
  brain[sp] <- list(neocortex$brain_mass_g_log_mean[i])
}

neocortex <- neocortex[match(artiodactyla$tip.label, neocortex$species), ]

m0 <- slouch.fit(
  phy,
  species = neocortex$species,
  response = neocortex$brain_mass_g_log_mean,
  fixed.fact = neocortex$diet,
  a_values = 0.01,
  sigma2_y_values = 1.0,
  anc_maps = "simmap",
  hillclimb = FALSE
)


theta <- m0$beta_primary$coefficients[,1]
lnl_brain_pruning <- sd_logL_pruning(phy, brain,
                                     sigma_squared = 1.0,
                                     alpha = 0.01, theta = theta)
lnl_brain_vcv <- m0$modfit$Support

test_that(paste("Simmap likelihood check"), {
  expect_equal(lnl_brain_pruning, lnl_brain_vcv)
})

