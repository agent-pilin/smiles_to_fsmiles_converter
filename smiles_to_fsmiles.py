# load library for regex
import re
import argparse

# this function takes as input sdf file and returns a string of chemical formula of the molecule in SMILES notation
def find_smiles_in_sdf(sdf_file):
    
    # open sdf file
    with open(sdf_file, 'r') as file:
        file_content = file.read()
    
    # define pattern corresponding to SMILES tag in sdf file
    smiles_tag = r'<PUBCHEM_OPENEYE_CAN_SMILES>\n(.*?)\n\n'

    # search for the SMILES notation formula with corresponding tag
    smiles_match = re.search(smiles_tag, file_content)

    # check if the SMILES tag is found
    if smiles_match:
        smiles_string = smiles_match.group(1)
        return smiles_string
    else:
        print("SMILES tag not found in the file.")
        return None

# this function takes as input SMILES string and returns list of SMILES substrings corresponding to chemical groups (rings or acyclic fragments)
def cut_smiles_into_groups(smiles_string):
    
    list_groups = [] # empty list
    ring_group_index = 1 # ring group index of (in SMILES group index starts with '1')
    
    # cut groups from SMILES strings while end of the string is not reached
    while len(smiles_string) > 0:
        # search for ring group pattern at the beginning of SMILES string ('C1', 'C2', 'N2', etc.)
        ring_group_match = re.search(r'^[A-Z][a-z]?{}'.format(ring_group_index), smiles_string)
        
        # if current group is ring
        if ring_group_match:
            # get current ring group
            group = re.search(r'^[A-Z][a-z]?{}.*'.format(ring_group_index) + r'[A-Z][a-z]?{}\)?'.format(ring_group_index), smiles_string).group()
            # remove current ring group from SMILES string
            smiles_string = re.sub(r'^[A-Z][a-z]?{}.*'.format(ring_group_index) + r'[A-Z][a-z]?{}\)?'.format(ring_group_index), '', smiles_string)
            # increment group index
            ring_group_index += 1
        # if current group is acyclic fragment
        else:
            # get current acyclic fragment
            group = re.sub(r'(^.*?)([A-Z][a-z]?{}.*)'.format(ring_group_index), r'\1', smiles_string)
            # remove current acyclic fragment from SMILES string
            smiles_string = re.sub(r'^.*?([A-Z][a-z]?{}.*)'.format(ring_group_index), r'\1', smiles_string)

        # add current group to the list
        list_groups.append(group)

        # if current group is the last one, stop while loop
        if group == smiles_string:
            break
        
    return list_groups

# this function takes as input chemical group and return its ring size ('0' for acyclic fragments)
def find_ring_size(group):
    
    # get ring group index from the group
    group_index = re.sub(r'^[A-Z][a-z]?([0-9]+).*', r'\1', group)

    # if current group is ring, ring group index will be found at least (?) twice
    if len(re.findall(r'[A-Z][a-z]?{}'.format(group_index), group)) > 1:
        # get ring group = substring from ring group starting with initial element (having group index) and ending with last element (having the same group index)
        ring_group = re.search(r'^[A-Z][a-z]?([0-9]+).*[A-Z][a-z]?\1', group).group()
        # remove subgroups from the ring group
        ring_group = re.sub(r'\(.*?\)', '', ring_group)
        # ring size is element count in ring group
        ring_size = len(re.findall('[A-Z][a-z]?', ring_group))
    # else current group is acyclic fragment
    else:
        ring_size = 0

    return ring_size

# check if chemical group is saturated or aromatic (all 'C' elements has one '=' bond)
def is_saturated(group):
    saturated = re.search(r'[^=\(\)][\(\)]?C[0-9]*[\(\)]?[^\(\)=0-9]', group)
    if saturated:
        return True
    else:
        return False
    
# this function takes as input list of chemical groups of the molecule in SMILES notation
# and returns the index of acyclic group within molecule (index in the list of chemical groups, NOT index corresponding to SMILES notation)
def find_acyclic_group(group_list):
    
    # by default index = 0, meaning that acyclic group within molecule not found
    acyclic_group_found = 0

    # iterate by group indices and groups
    for group_index, group in enumerate(group_list):
        # get group ring size
        ring_size = find_ring_size(group)
        # if acyclic group still not found
        if (acyclic_group_found == 0):
            # if current group is not marginal and acyclic
            if (ring_size == 0) & (group_index) & (group_index != len(group_list) - 1):
                # update the index of acyclic group
                acyclic_group_found = group_index
        # if acyclic group found, stop while loop
        else:
            break
            
    return acyclic_group_found


# this function takes as input chemical group in SMILES notation and returns itself but with element order counting in the opposite way (f.e. clockwise instead of counterclockwise)
def reverse_group_indexation(group):
    
    # get subgroup of chemical group if exists
    subgroup_size = 0
    subgroup_match = re.search(r'\(.*\)', group)
    if subgroup_match:
        subgroup = subgroup_match.group()
        # count element number in the subgroup
        subgroup_size = len(re.findall(r'[A-Z][a-z]?', subgroup))
    # count element number in the chemical group
    group_size = len(re.findall(r'[A-Z][a-z]?', group))

    # remove parentheses
    group = re.sub(r'[\(\)]', '', group)
    
    # cut chemical gorup into list of elements
    element_list = []
    while len(group) > 0:
        element_match = re.search(r'^[A-Z=][a-z]?[0-9]*', group)
        if element_match:
            element = element_match.group()
            # add element to the element list
            element_list.append(element)
            # remove element from the chain
            group = re.sub('^[A-Z=][a-z]?[0-9]*', '', group)

    group_rev = ''.join([element for element in element_list[::-1]])
    group_rev = re.sub(r'([A-Z][a-z]?[0-9]*)(=)([A-Z][a-z]?[0-9]*)', r'\1\3\2', group_rev)
    group_rev = re.sub(r'((=?[A-Z][a-z]?[0-9]*=?){' + str(group_size - subgroup_size - 2) + '}$)', r'(\1)', group_rev)
    return group_rev


# this function takes input a part of group (main of subgroup) and returns itself reversed
def reverse_chain(chain, group_index, chain_subgroup_match, subgroup):
    
    element_list = []
    # get elements from the chain one by one and add them to element list
    while len(chain) > 0:
        element_match = re.search(r'^[A-Z=][a-z]?[0-9]*', chain)
        if element_match:
            element = element_match.group()
            # remove group index from element if exists
            element = re.sub('[0-9]*', '', element)
            # add element to the element list
            element_list.append(element)
            # remove element from the chain
            chain = re.sub('^[A-Z=][a-z]?[0-9]*', '', chain)
    # reverse order of elements in the chain
    chain_rev = ''.join([element for element in element_list[::-1]])
    # if current chain is main chain
    if (subgroup == False):
        # add group index to the first element of the chain
        chain_rev = re.sub(r'(^[A-Z][a-z]?)', r'\g<1>{}'.format(group_index), chain_rev)
    # if chemical group does not have subgroup or current chain is subgroup
    if (chain_subgroup_match == False) | (subgroup == True):
        # add group index to the last element of the chain
        chain_rev = re.sub(r'([A-Z][a-z]?)(=?)\)?$', r'\g<1>{}\2'.format(group_index), chain_rev)
        
    return chain_rev

# this function takes as input SMILES formula of one chemical group and return its SMILES formula so that initial element is at the other connection point of chemical group)
def reverse_group(group):
    # get ring group index
    group_index = re.sub(r'^[A-Z][a-z]?([0-9]*).*', r'\1', group)
    # check if chemical group has a subgroup
    chain_subgroup_match = re.search(r'\(.*\)', group)
    if chain_subgroup_match:
        # get subgroup
        chain_subgroup = chain_subgroup_match.group()
        # remove paretheses
        chain_subgroup = re.sub('^\(', '', chain_subgroup)
        chain_subgroup = re.sub('\)$', '', chain_subgroup)
        # get main chain
        group = re.sub(r'\(.*\)', '', group)
        # reverse subgroup chain
        chain_subgroup_rev = reverse_chain(chain_subgroup, group_index, chain_subgroup_match, subgroup = True)
    # reverse main chain
    chain_main_rev = reverse_chain(group, group_index, chain_subgroup_match, subgroup = False)
    # define reversed group as main chain
    group_rev = chain_main_rev
    # add subgroup chain if exists
    if chain_subgroup_match:
        group_rev += '(' + chain_subgroup_rev + ')'
    
    return group_rev


# this function takes as input list of chemical groups and returns rearranged list so that acyclic group was in the first position
def rearrange_group_list(group_list, acyclic_group_found):
    
    # if acyclic group found
    if acyclic_group_found:
        # empty list that will be filled with group on one side of acyclic group
        group_list_rev = []
        # group indexing starts with the longest chain of groups counting from acyclic group
        if (acyclic_group_found > len(group_list) - 1):
            for group in group_list[(acyclic_group_found )::-1]:
                # reverse indexation for ring groups
                ring_size = find_ring_size(group)
                if ring_size:
                    group = reverse_group_indexation(group)
                # reverse group (like initial element is on the other side of the chain)
                group_rev = reverse_group(group)
                # add reversed group to the list
                group_list_rev.append(group_rev)
            # rewrite chain of chemical group on one side of acyclic group in initial list
            group_list[:(acyclic_group_found + 1):] = group_list_rev
        else:
            for group in group_list[(acyclic_group_found - 1)::-1]:
                # reverse indexation for ring groups
                ring_size = find_ring_size(group)
                if ring_size:
                    group = reverse_group_indexation(group)
                # reverse group (like initial element is on the other side of the chain)
                group_rev = reverse_group(group)
                # add reversed group to the list
                group_list_rev.append(group_rev)
            # rewrite chain of chemical group on one side of acyclic group in initial list
            group_list = group_list[(acyclic_group_found )::] + group_list_rev
            
    return group_list

# this function takes as input chemical group and its ring size in SMILES notation and returns its FSMILES notation
def transform_group(group, ring_size):

    # check if chemical group has multiple subgroups
    multiple_subgroups = len(re.findall('\)\(', group))
    if multiple_subgroups:
        # get first subgroup
        first_subgroup = re.search(r'\(.*?\)', group).group()
        # reverse first subgroup
        rev_first_subgroup = ''.join([x for x in first_subgroup[-2:0:-1]])
        # switch base element and first su group
        group = re.sub(r'([A-Z][a-z]?)\(.*?\)', r'{}\1()'.format(rev_first_subgroup), group)
        
    saturated = is_saturated(group)

    # remove double bond between C and N elements
    group = re.sub(r'([CN])([0-9]*)=([CN])', r'\1\2\3', group)
    # remove '=' at the beginning of subgroup
    group = re.sub(r'\(=C', r'(C', group)
    # add ring size for each element
    group = re.sub(r'=([A-Z][a-z]?)', r'=_{}\1'.format(ring_size), group)
    # add [*] at the end of the 
    group = re.sub('([A-Z][a-z]?[0-9]+$)', r'\1[*]', group)
    # remove ')' and '=)'
    group = re.sub('=?\)', '', group)
    # replace the start of subgroup by '([*])'
    group = re.sub('\(', '([*])', group)
    # add '_0' to initial element
    group = re.sub(r'([A-Z][a-z]?)[0-9]+', r'\g<1>1_0', group)
    # add suffix with ring size
    group = re.sub(r'([A-Z])', r'\1_{}'.format(ring_size), group)
    # add '_0' to '([*])'
    group = re.sub(r'(\(?\[\*\]\)?)', r'\1_0', group)

    # if group is aromatic or acyclic, then lowercase their elements
    if (saturated == False) & (ring_size > 0):
        group = group.lower()

    return group


# this function takes as input the list of chemical groups in FSMILES notation and returns molecule formula in FSMILES notation
def concat_groups(group_list):
    
    # define FSMILES string with starting tag
    fsmiles_string = "'start_0'"

    # iterate by chemical groups
    for group in group_list:
        # elongate FSMILES string with current group and separator tag
        fsmiles_string += group + "'sep_0’"

    # add ending tag
    fsmiles_string += "'end_0'"
    
    return fsmiles_string

# main function
def main():
    
    # create ArgumentParser object
    parser = argparse.ArgumentParser(description = 'SMILES to FSMILES conversion')

    # add arguments
    parser.add_argument('sdf_filepath', type = str, help = 'sdf filepath')
    parser.add_argument('fsmiles_filepath', type = str, help = 'filepath to write fsmiles string')

    # parse the command line arguments
    args = parser.parse_args()
    
    # declare arguments
    sdf_filepath = args.sdf_filepath
    fsmiles_filepath = args.fsmiles_filepath

    # get SMILES string from sdf file
    smiles_string = find_smiles_in_sdf(sdf_filepath)

    # cut SMILES into list of chemical groups
    group_list = cut_smiles_into_groups(smiles_string)

    # find non-marginal acyclic group
    acyclic_group_found = find_acyclic_group(group_list)

    # rearrange list of chemical groups
    group_list = rearrange_group_list(group_list, acyclic_group_found)

    # SMILES to FSMILES transformation by chemical group
    list_transformed = []
    for group in group_list:
        ring_size = find_ring_size(group)
        saturated = is_saturated(group)
        transformed_group = transform_group(group, ring_size)
        list_transformed.append(transformed_group)

    # concatenate FSMILES chemical groups
    fsmiles_string = concat_groups(list_transformed)

    # write FSMILES into text file
    with open(fsmiles_filepath, 'w') as file:
        file.write(fsmiles_string)

if __name__ == "__main__":
    main()
