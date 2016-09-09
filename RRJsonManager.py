
# import classes
import json

class RRJsonManager:
 
  # constructor
  def __init__(self):
    pass
  
  # read contacts in json format 
  def readConts(self, jsonfile):
    ''' jsonfile - file in jason format with data 
	containing contacts either in protein target or in prediction
	the data structures for these two cases are different
    '''
    with open(jsonfile) as f:
      data = json.load(f)
      f.close()
      return data
  
  # reformat targets contacts elemenibating the residues names
  def reformatTargetConts(self, data):
    result = {}
    for r1 in data:
      #r1_ = r1[3:] # obsolete format
      r1_ = r1
      lr2_ = []
      for r2 in data[r1]:
        #r2_ = r2[3:] # osolete format
        r2_ = r2
        lr2_.append(r2_)
      result[r1_] = lr2_
    return result
  
  # filter prediction contacts by probabilty
  def filterRRcontacts(self, rr_data, targetConts, prob=0.0, rr_len = "FL"):
    if rr_len not in ("10", "L5", "L2", "FL"):
      rr_len = "FL"
    result = []
    for chain in rr_data:
      for cont in rr_data[chain][rr_len]:
        if float(cont[2]) >= prob:
          try:
            # check and adjust enumeration
            if cont[1] in targetConts[cont[0]]:
              result.append( (min(int(cont[0])-1, int(cont[1])-1), max(int(cont[0])-1, int(cont[1]))-1) )
          except KeyError:
            pass
      break
    return result


# for testing/debigging
if __name__ == "__main__":
   targetFile = "input/T0859/T0859_contactsML.json";
   modelFile = "input/T0859/T0859RR287_1_trimmedLists_ML.json";
   jm = RRJsonManager()
   targetConts = jm.reformatTargetConts(jm.readConts(targetFile))
   rr_data = jm.readConts(modelFile)
   rr_conts = jm.filterRRcontacts(rr_data, targetConts, rr_len = "L5")
   print len(rr_conts)
   for cont in rr_conts:
     print cont[0], cont[1]
     break

