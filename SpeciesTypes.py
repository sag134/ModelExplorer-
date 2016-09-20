import sys
import xml.etree.ElementTree as ET
from converter import toString

def SpeciesTypes(tree):
	root = tree.getroot()
	el = [x for x in root[0]]
	elTag = [x.tag for x in root[0]]
	str1 = 'listOfSpeciesTypes'
	match = [i for i in range(0,len(elTag)) if str1 in elTag[i]]
	if len(match)>1:
		print 'Problem. >1 match. Breaking'
		sys.exit()
	else:
		match = match[0]
	speciesTypes = el[match]
	str2 = 'speciesType'
	str3 = 'bindingSiteSpeciesType'
	species = [x for x in speciesTypes if str2 in x.tag]
	bindingsites = [x for x in speciesTypes if str3 in x.tag]
	prefix = '{http://www.sbml.org/sbml/level3/version1/multi/version1}' #Do something about this. Instead of hard coding it.
	bindingsitenames = [dict(a.attrib)[prefix+'name'] for a in bindingsites]
	bindingsiteID = [dict(a.attrib)[prefix+'id'] for a in bindingsites]
	#print bindingsitenames
	#print bindingsiteID
	counter = 0
	molecules = [x for x in species]# if all(['listOfInSpeciesTypeBonds' not in i.tag for i in x])]
	sp = {}
	molID ={}
	cmplx = {}
	id2name = {}
	id2complexname = {}
	Molecules = {}
	for i in molecules:
		if all(['listOfInSpeciesTypeBonds' not in x.tag for x in i]):
			molecule_id = dict(i.attrib)[prefix+'id']
			if dict(i.attrib)[prefix+'name']:
				molecule_name = dict(i.attrib)[prefix+'name']
			else:
				molecule_name = ''
			full_molecule = (molecule_id,molecule_name)
			Molecules[full_molecule] = {'SpeciesTypes':[],'FeatureTypes':{}}
			for j in i:
				if prefix+'listOfSpeciesTypeInstances' in j.tag:
					for k in j:
						id_ = dict(k.attrib)[prefix+'speciesType']
						if dict(k.attrib)[prefix+'name']:
							name_ = dict(k.attrib)[prefix+'name']
						else:
							name_ = ''
						Molecules[full_molecule]['SpeciesTypes'].append([id_,name_])
						#print id_
				if prefix+'listOfSpeciesFeatureTypes' in j.tag:
					for k in j:
						id_ = dict(k.attrib)[prefix+'id']
						if dict(k.attrib)[prefix+'name']:
							name_ = dict(k.attrib)[prefix+'name']
						else:
							name_ = ''
						feature = (id_,name_)
						Molecules[full_molecule]['FeatureTypes'][feature] = []
						for q in k:
							if 'listOfPossibleSpeciesFeatureValues' in q.tag:
								for q1 in q:
									id_ = dict(q1.attrib)[prefix+'id']
									if dict(q1.attrib)[prefix+'name']:
										name_ = dict(q1.attrib)[prefix+'name']
									else:
										name_ = ''	
									Molecules[full_molecule]['FeatureTypes'][feature].append([id_,name_])					
	nmol = len(Molecules)
	keys = Molecules.keys()
	for i in range(0,nmol):
		key = keys[i]
		st = Molecules[key]['SpeciesTypes']
		nst = len(st)
		for j in range(0,nst):
			spt = st[j]
			#Check if the species type is a binding site
			spt_id = spt[0]
			if spt_id not in bindingsiteID:
				#spt_id must appear as a key 
				tkey = [x[0] for x in keys]
				if spt_id not in tkey:
					print 'Something is wrong'
					return
				index = tkey.index(spt_id)
				#Replace j in st with species types
				stnew = Molecules[keys[index]]['SpeciesTypes']
				Molecules[key]['SpeciesTypes'][j] = []
				for x in stnew:
					Molecules[key]['SpeciesTypes'].append(x)
				ftnew = Molecules[keys[index]]['FeatureTypes']
				Molecules[key]['FeatureTypes'].update(ftnew)

		Molecules[key]['SpeciesTypes'] = [x for x in Molecules[key]['SpeciesTypes'] if len(x)>0]	
	Complexes = {}
	for i in molecules:
		if not(all(['listOfInSpeciesTypeBonds' not in x.tag for x in i])):
			complex_id = dict(i.attrib)[prefix+'id']
			if dict(i.attrib)[prefix+'name']:
				complex_name = dict(i.attrib)[prefix+'name']
			else:
				complex_name = ''
			full_complex = (complex_id,complex_name)
			Complexes[full_complex] = {'molecules':[],'bonds':[],'components':[]}
			for j in i:
				if prefix+'listOfSpeciesTypeInstances' in j.tag:
					for k in j:
						id_ = dict(k.attrib)[prefix+'speciesType']
						if dict(k.attrib)[prefix+'name']:
							name_ = dict(k.attrib)[prefix+'name']
						else:
							name_ = ''		
						Complexes[full_complex]['molecules'].append([id_ ,name_])
				if prefix+'listOfInSpeciesTypeBonds' in j.tag:
					for k in j:
						bs1=dict(k.attrib)[prefix+'bindingSite1']
						bs2=dict(k.attrib)[prefix+'bindingSite2']
						Complexes[full_complex]['bonds'].append([bs1,bs2])
				if prefix+'listOfSpeciesTypeComponentIndexes' in j.tag:
					for k in j:
						id_ =dict(k.attrib)[prefix+'id']
						comp_=dict(k.attrib)[prefix+'component']
						Complexes[full_complex]['components'].append([id_,comp_])

			#Check if the bond definition is in the components list. If so, replace the bond id, with the component name.
			for component in Complexes[full_complex]['components']:
				component_name = component[0]
				bondcounter = 0
				for bond in Complexes[full_complex]['bonds']:
					if component_name == bond[0]:
						Complexes[full_complex]['bonds'][bondcounter][0] = component[1]
					elif component_name == bond[1]:
						Complexes[full_complex]['bonds'][bondcounter][1] = component[1]
					bondcounter = bondcounter + 1


	str1 = 'listOfSpecies'
	prefix1 ='{http://www.sbml.org/sbml/level3/version1/core}'
	match = [i for i in range(0,len(elTag)) if prefix1+str1 == elTag[i]]
	if len(match)>1:
		print 'Problem. >1 match. Breaking'
		sys.exit()
	else:
		match = match[0]
	speciesList = el[match]
	ModelSpecies = {}
	for sp in speciesList:
		if 'species' in sp.tag:
			id_ = dict(sp.attrib)['id']
			ModelSpecies[id_] = {'name':dict(sp.attrib)['name'],'speciesType':dict(sp.attrib)[prefix+'speciesType'],'bindingsites':[],'features':[]}
			if 'initialConcentration' in dict(sp.attrib).keys():
				ModelSpecies[id_]['ic']=dict(sp.attrib)['initialConcentration']
		for j in sp:
			if 'BindingSites' in j.tag:
				for q in j:
					if dict(q.attrib)[prefix+'bindingStatus'] != 'either':
						ModelSpecies[id_]['bindingsites'].append(dict(q.attrib)[prefix+'component'])
			if 'SpeciesFeatures' in j.tag:
				for q in j:
					ModelSpecies[id_]['features'].append([dict(q.attrib)[prefix+'speciesFeatureType'],dict(q[0][0].attrib)[prefix+'value']])
	Result = {'Molecules':Molecules,'Complexes':Complexes,'Species':ModelSpecies}
	return Result


def FlattenDict(x):
	for i in x.keys():
		if x[i] in x.keys():
			x[i] = x[x[i]]
	return x			

if __name__ == '__main__':
	#TESTING
	#Parse sbml file
	tree = ET.parse('../Downloads/simple_system_sbml_sbmlmulti.xml')
	#tree = ET.parse('ex15.xml')
	moleculetypes = SpeciesTypes(tree)