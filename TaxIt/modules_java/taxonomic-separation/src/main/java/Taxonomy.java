/* Copyright (c) 2016,
 * Mathias Kuhring, KuhringM@rki.de, Robert Koch Institute, Germany,
 * All rights reserved. For details, please note the license.txt.
 */

import java.util.*;
import java.util.Map.Entry;

/**
 * Class for the integration of all taxonomic information (ID mapping, NCBI taxonomix nodes and names)
 *
 */
public class Taxonomy {

	Nodes nodes;
	Names names;
	IDMapping idmapping;

	// the "active" nodes (taxid with corresponding gis) have 2 additional Maps:
	// taxidToLineage contains the path of node to the root as a list of taxids.
	// each list starts with the root and ends with the node in question.
	private HashMap<Integer,ArrayList<Integer>> taxidToLineage;
	// nodeToDescendants contains a set of all descendants of node (no particular order).
	// it represents a subtree, including the root node as well.
	private HashMap<Integer,HashSet<Integer>> nodeToDescendants;


	public Taxonomy(String idMappingFile,
					String taxNodesFile,
					String taxNamesFile){

		idmapping = new IDMapping();
		idmapping.importMapping(idMappingFile);
		int totalTax = idmapping.countTaxids();
		System.out.println(totalTax + " orgs in total... ");

		System.out.println("tax nodes and names import... ");
		nodes = new Nodes();
		nodes.parseNodes(taxNodesFile);

		idmapping.removeTaxidsWithoutNodes(
				new HashSet<Integer>(nodes.getTaxidToParent().keySet()));

		names = new Names();
		names.parseNames(taxNamesFile, idmapping.getAllTaxids());

		taxidToLineage = new HashMap<Integer,ArrayList<Integer>>();
		nodeToDescendants = new HashMap<Integer,HashSet<Integer>>();
	}

	public void calcLineage(HashMap<Integer, HashSet<String>> taxidToIdMap,
							Nodes nodes){
		HashMap<Integer, Integer> taxidToParent = nodes.getTaxidToParent();
		HashMap<Integer, String> taxidToRank = nodes.getTaxidToRank();
		HashMap<String, HashSet<Integer>> rankToTaxids = nodes.getRankToTaxids();

		HashSet<Integer> uniqueTaxids = new HashSet<Integer>(taxidToIdMap.keySet());

		ArrayList<Integer> newLineage;
		int lastTaxid, nextTaxid;
		String rank;

		for (Integer taxid : uniqueTaxids){
			if (!taxidToLineage.containsKey(taxid)){
				newLineage = new ArrayList<Integer>();
				newLineage.add(taxid);
				lastTaxid = taxid;

				// put own taxid in nodeToDescendants, necessary e.g. when sampling from a branch
				if (!this.nodeToDescendants.containsKey(taxid)){
					this.nodeToDescendants.put(taxid, new HashSet<Integer>());
				}
				this.nodeToDescendants.get(taxid).add(taxid);

				while (lastTaxid != 1){
					nextTaxid = taxidToParent.get(lastTaxid);
					newLineage.add(0, nextTaxid);

					if (!this.nodeToDescendants.containsKey(nextTaxid)){
						this.nodeToDescendants.put(nextTaxid, new HashSet<Integer>());	
						this.nodeToDescendants.get(nextTaxid).add(nextTaxid);
					}
					this.nodeToDescendants.get(nextTaxid).addAll(nodeToDescendants.get(lastTaxid));

					rank = taxidToRank.get(lastTaxid);
					if (!rankToTaxids.containsKey(rank)){
						rankToTaxids.put(rank, new HashSet<Integer>());
					}
					rankToTaxids.get(rank).add(lastTaxid);

					lastTaxid = nextTaxid;
				}

				taxidToLineage.put(taxid, newLineage);
			}
		}
	}

	// assign proteins of a non-leaf node to all of its leafs
	public void reassignProteins(IDMapping mapping){
		HashMap<Integer, HashSet<String>> taxidToGi = mapping.getTaxidToId();
		HashMap<String, HashSet<Integer>> giToTaxid = mapping.getIdToTaxid();

		HashSet<Integer> subnodes;
		int taxid;
		Iterator<Entry<Integer, HashSet<String>>> it = taxidToGi.entrySet().iterator();
		while (it.hasNext()) {
			taxid = it.next().getKey();
			subnodes = nodeToDescendants.get(taxid);
			// if node has children
			if (subnodes.size() > 1){
				for (int subnode : subnodes){
					// assign proteins to children, if they are a leaf
					if (nodeToDescendants.get(subnode).size() == 1){
						taxidToGi.get(subnode).addAll(taxidToGi.get(taxid));
						for (String gi : taxidToGi.get(taxid)){
							giToTaxid.get(gi).add(subnode);
						}
					}
				}
				for (String gi : taxidToGi.get(taxid)){
					giToTaxid.get(gi).remove(taxid);
				}
				it.remove();
			}
		}
	}

	public void cleanIdMap(HashSet<String> ids) {
		System.out.println("clean up ids... ");
		ArrayList<String> removed = idmapping.removeIdsWithNoTaxidMapping(ids);
		if (removed.size() > 0){
			System.out.println("removed " + removed.size() + " ids without taxid mapping... ");
			System.out.println("\tfirst 3 removed ids: " + removed.subList(0, Math.min(3, removed.size())));
		}
	}

	public void reassignIds(int iteration) {
		System.out.println("calculate lineages... ");
		calcLineage(idmapping.getTaxidToId(), nodes);

		System.out.println("down propagation of redundant proteins... ");
		reassignProteins(idmapping);

		if (iteration == 2){
			System.out.println("up propagation of strain proteins... ");
			reassignStrainProteins(idmapping);
		}
	}

	public HashSet<Integer> getTaxids(HashSet<String> matches) {
		return idmapping.getTaxids(matches);
	}

	public HashMap<String, HashSet<Integer>> getIdTaxidMap(HashSet<Integer> taxids) {

		HashMap<String, HashSet<Integer>> idTaxidMap = new HashMap<String, HashSet<Integer>>();

		for (int taxid : taxids){
			if (idmapping.hasGis(taxid)){
				for (String gi : idmapping.getIds(taxid)){
					if (!idTaxidMap.containsKey(gi)){
						idTaxidMap.put(gi, new HashSet<Integer>());
					}
					idTaxidMap.get(gi).add(taxid);
				}
			}
		}

		return idTaxidMap;
	}

	// assign proteins of a non-leaf node to all of its leafs
	public void reassignStrainProteins(IDMapping mapping){
		HashMap<Integer, HashSet<String>> taxidToGi = mapping.getTaxidToId();
		HashMap<String, HashSet<Integer>> giToTaxid = mapping.getIdToTaxid();

		int pTaxid;
		HashSet<Integer> taxidsToRemove = new HashSet<Integer>();
		HashMap<Integer, HashSet<String>> taxidToGoToAdd = new HashMap<Integer, HashSet<String>>();

		for (int taxid : taxidToGi.keySet()){
			// if node has "species" parent
			if ((pTaxid = getRankParent(taxid, "species")) > 0){

				// assign proteins to "species" parent
				if (!taxidToGoToAdd.containsKey(pTaxid)){
					taxidToGoToAdd.put(pTaxid, new HashSet<String>());
				}
				taxidToGoToAdd.get(pTaxid).addAll(taxidToGi.get(taxid));

				for (String gi : taxidToGi.get(taxid)){
					giToTaxid.get(gi).add(pTaxid);
				}

				// remove protein taxid mappings
				for (String gi : taxidToGi.get(taxid)){
					giToTaxid.get(gi).remove(taxid);
				}
				// remember to remove taxid
				taxidsToRemove.add(taxid);
			}
		}

		// remove taxids lower than species
		taxidToGi.keySet().removeAll(taxidsToRemove);
		// add new species mappings
		taxidToGi.putAll(taxidToGoToAdd);
	}

	// check whether a node has a specific ranked parent and return the taxid, else 0
	private int getRankParent(int taxid, String rank){
		int pTaxid = 0;
		ArrayList<Integer> lineage = taxidToLineage.get(taxid);
		for (int i=0; i<lineage.size()-1; i++){ // skip last element, which is the node itself
			if (nodes.getRank(lineage.get(i)).matches(rank)){
				pTaxid = lineage.get(i);
				break;
			}
		}
		return pTaxid;
	}
}
