import numpy as np
import matplotlib.pyplot as plt

"""
A program to perform clustering of elements (i.e. site states) using the average absolute density matrix (adm).

The criterium for clustering two sites is a threshold value of their mutual off-diagonal adm element, containing information on
how strongly sites 'talk to each other', absorbing both their coupling and respective energy gap.
Two clusters (of sites) are coupled if at least one site in one cluster is 'talking sufficiently strongly' to any site in the other cluster.

This algorithm calculates how the number of indentifiable clusters depends on the above threshold value, yielding
    1) a plot showing the number of clusters for the specified range of epsilon (threshold)
    2) a 'dendrogram' type of plot depicting a more detailed way of the clustering for different threshold values
    3) a file 'segments_found.dat', specifying to which cluster each site belongs, at the given desired number of clusters,
       which can be used as input for MC-FRET in NISE_2017.
       NB: the above plots might help in choosing a 'relevant' number of clusters.

The input for this program is an (average) absolute density matrix that can be obtained with the 'Analyse' technique within NISE_2017.

@author: Gijsbert ten Hoven
"""

##### initialising #####

adm =  np.loadtxt("AbsoluteDensityMatrix.dat")
nr_sites = adm.shape[0]

#epsilons = np.linspace(1,0.3,1000)
epsilons = np.linspace(0,1,10)

# epsilons array must decrease: proper ordering needed
epsilons = np.sort(epsilons)[::-1]

# will store a list of lists (that contain all site
# indices in a segment) for each value of epsilon.
segments_epsilon = [] 

# will store the number of segments found at specified epsilons
nr_of_segments = [] 

# set initial conditions: all sites are in their own segment
segments = [[i] for i in range(nr_sites)]

##### define subroutines for the main clustering loop #####

def merge_routine(segments, epsilon):
    # double loop over existing segments
    for seg_i in range(len(segments)):
        for seg_j in range(len(segments))[seg_i+1:]:
            if connection_check(seg_i, seg_j, segments, epsilon):
                updated_segments =  merge_two_segments(seg_i, seg_j, segments)
                return True, updated_segments
    return False, segments

def connection_check(seg_i, seg_j, segments, epsilon):
    # check if two segments are connected
    for site_m in segments[seg_i]:
        for site_n in segments[seg_j]:
            off_diag = adm[site_m, site_n]
            diag_m = adm[site_m,site_m]
            diag_n = adm[site_n,site_n]
            if off_diag > epsilon * diag_m or off_diag > epsilon * diag_n:
                # no need for further checking, one connection is sufficient
                return True
    return False

def merge_two_segments(seg_i, seg_j, segments):
    """ merge two segments into one, return updated segment list """
    seg_A = segments[seg_i].copy()
    seg_B = segments[seg_j].copy()
    merged = seg_A + seg_B
    segments.remove(seg_A)
    segments.remove(seg_B)
    segments.append(merged)
    return segments

def write_segment_file(segments):
    "Write the identified segments in a file, using the input format for the MC-FRET technique in NISE_2017"
    with open("segments_found.dat", "w") as file:
        file.write(str(nr_sites) + '\n')
        for site_k in range(nr_sites):
            for segment_i, segment in enumerate(segments):
                if site_k in segment:
                    file.write(str(segment_i) + ' ')
    file.close()

##### Main loop ####

for epsilon in epsilons: # skip the initial condition
    print("Finding segments for epsilon {:3f}".format(epsilon), end='\r')
    # check whether any two segments can be merged
    # continue to do the merging routine until no further merges possible at this value of epsilon
    merge_possible = True # dummy initial value

    while merge_possible == True:
        # check if new merges are possible, until no further merges can be done at this value of epsilon
        merge_possible, segments = merge_routine(segments, epsilon)
        
    # no further merges possible: segments determined at this value of epsilon
    segments_epsilon.append(segments.copy()) # 'copy' in order to avoid storing a pointer
    nr_of_segments.append(len(segments))

    # bonus: make a segment.dat file as input for MC-FRET in NISE_2017
    # e.g. for the case that two segments are found, to check if this yields the two separate tubes in the chlorosomes
    # It can happen that the desired number of clusters is not found, in which case a more detailed range of epsilons may help, depending on the system investigated
    desired_segment_number = 2
    if len(segments) == desired_segment_number:
        write_segment_file(segments)

##### plot the number of segments as function of the epsilon parameter #####

fig0 = plt.figure()
ax0 = fig0.subplots(1)
ax0.plot(epsilons, nr_of_segments)
ax0.set_title("Segments in system as function of $\epsilon$ parameter in density marix")
ax0.set_xlabel("$\epsilon$ parameter")
ax0.set_ylabel("nr. of segments")

##### make a dendrogram: show the tree structure of the clustering as function of the epsilon parameter #####

fig1 = plt.figure()
ax = fig1.subplots(1)
ax.set_ylabel("$\epsilon$ parameter")
ax.set_title("Splitting of system in clusters depending on $\epsilon$\n Total elements: " + str(nr_sites))
ax.set_xticks([])
ax.set_xlim(-0.1, 1.13)
unit_width = 1 / nr_sites

for epsilon_i, epsilon in reversed(tuple(enumerate(epsilons))):
    # connect the clusters to their sub_clusters
    # after that is done, space the dendrogram out evenly over the plot
    n_segments = nr_of_segments[epsilon_i]
    epsilons_plot = np.ones(n_segments) * epsilon
    horizontal_positions = []
    cluster_widths = []
        # make the connections in the tree
    if epsilon_i != len(epsilons) - 1:
        # horizontal
        prev_epsilon = epsilons[epsilon_i + 1]
        # vertical 
        for segment in segments_epsilon[epsilon_i]:
            cluster_width = unit_width * len(segment)
            cluster_widths.append(cluster_width)
            # find the parent segment to which it belongs
            for prev_seg_j, prev_segment in enumerate(segments_epsilon[epsilon_i + 1]):
                if set(segment) <= set(prev_segment):
                    x = cluster_width / 2 + parent_ref_x[prev_seg_j]
                    horizontal_positions.append(x)
                    parent_ref_x[prev_seg_j] += cluster_width
                    ax.plot([x, parent_positions[prev_seg_j]],[epsilon, prev_epsilon], color = 'black',alpha=0.3)
                    break

    else:
        # the first set to put in the figure
        ref_x = 0
        for segment in segments_epsilon[epsilon_i]:
            cluster_width = unit_width * len(segment)
            cluster_widths.append(cluster_width)
            horizontal_positions.append(cluster_width / 2 + ref_x)
            ref_x += cluster_width
    
    ax.scatter(horizontal_positions, epsilons_plot, s = 0.3, color = 'black')
    ax.text(1.04, epsilon, str(nr_of_segments[epsilon_i]), fontsize='small', verticalalignment='center')
    parent_positions = horizontal_positions.copy()
    parent_widths = cluster_widths.copy()
    parent_ref_x = np.array(parent_positions) - np.array(parent_widths) / 2

plt.show()

