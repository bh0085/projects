"""### funzioneSpazzino ###
usage: listOutput = funzioneSpazzino(list1,list2)

where: * list1 = list of the unwanted strings that should be cancelled from each item of the list2
	   * list2 = list of the strings wich should be epurated from EVERY string listed in list1

eg.

>> import funzioneSpazzino
>> a = ["ab","cd"]
>> b = ["abtestcd", "testabcd", "teabstcdabcd"]
>> c = funzioneSpazzino.funzioneSpazzino(a,b)
>> print c
["test", "test", "test"]
>>

Claudio Virili     claudio.virili@gmail.com
Emanuele Casali    emanuele.casali@gmail.com

22/01/2010
"""

import re

def funzioneSpazzino(nonVolute,listaFrasi):
	listaFrasiPulite = []
	for frase in listaFrasi:
		# E' necessario ciclare tutte le parole su una singola frase !
		for parola in nonVolute:
		    patternObj = re.compile(parola)
		    x = 1
	            while x == 1:
				b = re.search(patternObj, frase)
				if type(b) == type(None):
					x = 0 # Condizione di uscita
				else:
					frase = frase[0:b.start()]+frase[b.start()+len(parola):len(frase)]
		listaFrasiPulite.append(frase)	
	return listaFrasiPulite


			
