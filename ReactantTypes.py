import sys
import xml.etree.ElementTree as ET
from converter import toString
from SpeciesTypes import *
import re

def Reactions(tree):
	prefix = '{http://www.sbml.org/sbml/level3/version1/core}'
	root = tree.getroot()
	el = [x for x in root[0]]
	elTag = [x.tag for x in root[0]]
	match = [i for i in range(0,len(elTag)) if 'listOfReactions' in elTag[i]]
	if len(match)>1:
		print 'Problem. >1 match. Breaking'
		sys.exit()
	else:
		match = match[0]
	rxns = el[match]
	reactants = {}
	products = {}
	parameters = {}
	nrxns = len(rxns)	
	#print nrxns
	ReactionList = []
	for i in range(0,nrxns):
		tmp = {}
		rxn = rxns[i]
		reactant = []
		if 'name' in dict(rxn.attrib).keys():
			tmp['annotation'] = dict(rxn.attrib)['name']
		else:
			tmp['annotation'] = ''
		tmp['reversible'] = dict(rxn.attrib)['reversible']
		tmp['reactants'] = []
		tmp['products'] = []
		tmp['parameters'] = []
		for j in rxn:
			if j.tag == prefix+'listOfReactants':
				for k in j:
					tmp['reactants'].append(dict(k.attrib)['species'])
			elif j.tag == prefix+'listOfProducts':
				for k in j:
					tmp['products'].append(dict(k.attrib)['species'])
			elif j.tag == prefix+'kineticLaw':
				for k in j:
					if k.tag == prefix+'listOfLocalParameters':
						for q in k:
							tmp['parameters'].append({dict(q.attrib)['id']:dict(q.attrib)['value']})
		ReactionList.append(tmp)
	return ReactionList


def rxnstring(RL,tree):
	s = SpeciesTypes(tree)['Species']
	m = SpeciesTypes(tree)['Molecules']
	c = SpeciesTypes(tree)['Complexes']
	result = PartialMolecule(s,m,c)
	pm = result['pm']
	pc = result['pc']
	string_pm = toString(pm)
	#print pm
	string_pc = Complex2String(pc)
	reactantpatterns = []
	productpatterns = []
	reversible = []
	for i in RL:
		rpattern = []
		ppattern = []	
		if i['reversible'] != 'false':
			reversible.append(' <-> ')
		else:
			reversible.append(' -> ')
		for j in i['reactants']:
			if j in string_pm.keys():
				rpattern.append(string_pm[j])
			elif j in string_pc.keys():
				rpattern.append(string_pc[j])
		for j in i['products']:
			if j in string_pm.keys():
				ppattern.append(string_pm[j])
			elif j in string_pc.keys():
				ppattern.append(string_pc[j])
		if len(i['products']) == 0:
			ppattern = ['0']
		reactantpatterns.append(rpattern)
		productpatterns.append(ppattern)
	rxnstring = 'begin reaction rules \n'

	for i in range(0,len(RL)):
		if RL[i]['annotation'] == '':
			RL[i]['annotation'] = 'R'+str(i)
		rxnstring = rxnstring + RL[i]['annotation']+': '
		for pattern in reactantpatterns[i]:
			rxnstring = rxnstring+pattern+' + '
		rxnstring = rxnstring[:-2]+reversible[i]
		for pattern in productpatterns[i]:
			rxnstring = rxnstring+pattern+' + '
		rxnstring = rxnstring[:-2]+'\t'
		for parameters in RL[i]['parameters']:
			for key in parameters.keys():
				rxnstring = rxnstring+key+','
		rxnstring = rxnstring[:-1]+'\n'
	rxnstring = rxnstring+'end reaction rules \n'
	return rxnstring

def toString(Molecules):
	mol_string = {}
	for i in  Molecules:
		mol_name = i[1]
		m = mol_name+'('
		st = Molecules[i]['SpeciesTypes']
		ft = Molecules[i]['FeatureTypes']
		#print ft
		list_of_st = []
		flag = 0
		for j in st:
			if j[1] != '':
				m = m+j[1] +','
			else:
				m = m + j[0] +','
			flag = 1
		for j in ft.keys():
			m = m+j[1]
			if any(isinstance(el, list) for el in ft[j]):
				for k in ft[j]:
					m = m+'~'+k[1]
			else:
				m = m +'~'+ft[j][1]
			m = m+','
			flag = 1
		if flag == 0:
			m = m+')'
		else:
			m = m[:-1]+')'
		mol_string[i[0]] = m
	#print mol_string
	return mol_string

def Complex2String(complexes):
	complex_string = {}
	for i in complexes:
		c = ''
		for j in complexes[i]['molecules']:
			mname = toString(j).values()[0]
			mol = j.keys()[0][1]
			mol  = mol.replace('(','\(')
			mol = mol.replace(')','\)')
			counter = 0
			bondCounter = {}
			for b in complexes[i]['bonds']:
				counter = counter + 1
				st = j[j.keys()[0]]['SpeciesTypes']
				bondCounter[(b[0][0],b[1][0])] = counter
				for bond in b:
					if bond in st:
						bst = bond[1]
						m = re.search(mol+'\(\s*'+bst+'\s*[,\)]',mname)
						index = mname.find(m.group(0))+len(m.group(0))
						tmp = bondCounter.keys()
						ind = 1e29
						for tmp1 in tmp:
							if bond[0] == tmp1[0] or bond[0] == tmp1[1]:
								ind = bondCounter[tmp1]
								break
						mname_new = mname[0:index-1]+'!'+str(ind)+mname[index-1:]
			c = c+mname_new+'.'
		c = c[:-1]
		complex_string[i[0]] = c
	return complex_string

def PartialMolecule(Species,Molecules,Complexes):
	complex_id = [x[0] for x in Complexes.keys()]
	mol_id = [x[0] for x in Molecules.keys()]
	partialMolecule = {}
	partialComplex = {}
	for i in Species:
		st = Species[i]['speciesType']
		f = Species[i]['features']
		bs = Species[i]['bindingsites']
		pf_id = [x[0] for x in f]
		pf_val_id = [x[1] for x in f]
		if st in mol_id:
			mind = mol_id.index(st)
			key = Molecules.keys()[mind]
			pname = (i,key[1])
			partialMolecule[pname] = {'SpeciesTypes' : [], 'FeatureTypes':{}}
			full_molecule = Molecules[key]
			full_features = full_molecule['FeatureTypes']
			for ff in full_features.keys(): #Complete feature list
				fid = ff[0] #Check if feaure id is present
				if fid in pf_id:
					fval = full_features[ff] #Complete feature values list
					for j in fval:
						if j[0] in pf_val_id:
							partialMolecule[pname]['FeatureTypes'][ff] = j
			
			full_binding = full_molecule['SpeciesTypes']
			for binding in full_binding:
				if binding[0] in bs:
					partialMolecule[pname]['SpeciesTypes'].append(binding)	
		if st in complex_id:
			cid = complex_id.index(st)
			key = Complexes.keys()[cid]
			cname = (i,key[1])
			pmolecules = Complexes[key]['molecules']
			#Go through participating molecules
			features = Species[i]['features']
			ft = [x[0] for x in features]
			partialComplex[cname] = {'molecules':[],'bonds':[]}
			bst = Species[i]['bindingsites']
			bt = [x[0] for x in bst]
			bonds = Complexes[key]['bonds']
			for j in pmolecules:
				pid = mol_id.index(j[0])
				molecule =  Molecules[Molecules.keys()[pid]]
				moleculeCopy = {}
				moleculeCopy[Molecules.keys()[pid]] = {'SpeciesTypes':[],'FeatureTypes':{}}
				#Edit the molecule based on features provided
				mft =  [x[0] for x in molecule['FeatureTypes'].keys()]
				mst = [x[0] for x in molecule['SpeciesTypes']]
				for f in ft:
					if f in mft:
						fid = mft.index(f)
						sft = molecule['FeatureTypes'].keys()[fid]
						fvalues = molecule['FeatureTypes'][sft]
						fval = features[ft.index(f)][1]
						newfval = [x for x in fvalues if x[0] == fval]

						moleculeCopy[Molecules.keys()[pid]]['FeatureTypes'][sft] = newfval[0]
				for bond in bonds:
					b1 = bond[0]
					b2 = bond[1]
					if b1 in mst:
						b_index = mst.index(b1)
						site = molecule['SpeciesTypes'][b_index]
					if b2 in mst:
						b_index = mst.index(b2)
						site = molecule['SpeciesTypes'][b_index]
					if site not in moleculeCopy[Molecules.keys()[pid]]['SpeciesTypes']:
						moleculeCopy[Molecules.keys()[pid]]['SpeciesTypes'].append(site)
					#	partialComplex[cname]['bonds'].append(site)
				partialComplex[cname]['molecules'].append(moleculeCopy)
			for bond in bonds:
				b1 = bond[0]
				b2 = bond[1]
				tmp = []
				for j in pmolecules:
					pid = mol_id.index(j[0])
					molecule =  Molecules[Molecules.keys()[pid]]
					mst = [x[0] for x in molecule['SpeciesTypes']]
					if b1 in mst:
						b_index = mst.index(b1)
						site = molecule['SpeciesTypes'][b_index]
					if b2 in mst:
						b_index = mst.index(b2)
						site = molecule['SpeciesTypes'][b_index]
					tmp.append(site)
				partialComplex[cname]['bonds'].append(tmp)
	return {'pm':partialMolecule,'pc':partialComplex}

def MolTypesString(m):
	s = toString(m)
	moltypes = 'being molecule types \n'
	for j in s:
		moltypes = moltypes + s[j]+'\n'
	moltypes = moltypes + 'end molecule types \n'
	return moltypes

if __name__ == '__main__':
	#TESTING
	#Parse sbml file
	#tree = ET.parse('ecad.xml')
	#tree = ET.parse('ex15.xml')
	#tree = ET.parse('../Downloads/simple_system_sbml_sbmlmulti.xml')
	tree = ET.parse('/Users/sanjanagupta/Desktop/multi/bionetgen/bng2/toy-jim_sbml_sbmlmulti.xml')
	RL = Reactions(tree)
	s = SpeciesTypes(tree)['Species']
	m = SpeciesTypes(tree)['Molecules']
	c = SpeciesTypes(tree)['Complexes']
	#print toString(m)
	#print toString(m)
	x1 = MolTypesString(m)
	x2 = rxnstring(RL,tree)
	print x1
	print x2
	#result = PartialMolecule(s,m,c)
	#print result['pm'][result['pm'].keys()[1]]
	#print result['pc']
	#molstring = toString(result['pm'])
	#print molstring
	#cstring = Complex2String(result['pc'])
	#print cstring
	'''for i in RL:
		for j in i['reactants']:
			if j in molstring.keys():
				print molstring[j]'''
			#Construct partial molecule type
#			print s[j]
#		print i['reactants']

	#moleculetypes =  MoleculeTypes(tree)['cmplx']
	#for i in moleculetypes.keys():
	#	print toString(moleculetypes[i])
	#print RL[0]['reactants'][0]
	#rstring = rxnstring(RL)
	#print rstring



