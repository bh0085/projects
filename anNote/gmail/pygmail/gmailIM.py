""" ### gmail web engine ###

Claudio Virili     claudio.virili@gmail.com
Emanuele Casali    emanuele.casali@gmail.com

22/01/2010

"""


import mechanize
import cookielib
from BeautifulSoup import BeautifulSoup
import html2text
import sys
import os
import re
import funzioneSpazzino

sys.stdout.write(os.popen('clear').read()) # Pulisce lo schermo

	
# Inizio funzione per la "ricezione" delle mail
class retriveMail:
	def __init__(self, email, passwd):
		# Apertura del browser
		self.br = mechanize.Browser()
		self.cj = cookielib.LWPCookieJar()
		self.br.set_cookiejar(self.cj)
		self.br.set_handle_equiv(True)
		self.br.set_handle_redirect(True)
		self.br.set_handle_referer(True)
		self.br.set_handle_robots(False)
		self.br.set_handle_refresh(mechanize._http.HTTPRefreshProcessor(), max_time=1)
		self.br.addheaders = [('User-agent', 'Firefox/3.0.1')]

		# Login in gmail.com
		self.br.open('http://gmail.com')
		self.br.select_form(nr=0)
		self.br.form['Email'] = email
		self.br.form['Passwd'] = passwd
		sys.stdout.write(os.popen('clear').read()) # Pulisce lo schermo
		self.br.submit()
	
		# Cerca il link alla versione HTML di base e lo segue.
		linkHtmlbase = self.br.links(url_regex = '\?ui=html&zy=e')
		linkHtml = [y for y in linkHtmlbase]
		if linkHtml != []:
			self.br.follow_link(linkHtml[0])
		
	
	def retriveMainInfos(self, parametroPagina):
	# Creazione dell'obj soup dalla risposta della prima pagina (50 mail)
		zuppaPrimaPagina = BeautifulSoup(self.br.response().read())
		if parametroPagina == "1": 
			self.br.open(self.funzLinkBase(zuppaPrimaPagina)+"?s=i") # Se e' 1 va all'inbox
		elif parametroPagina == "2":
			self.br.open(self.funzLinkBase(zuppaPrimaPagina)+"?s=s") # Se e' 2 va al sent
		elif parametroPagina == "3":
			self.br.open(self.funzLinkBase(zuppaPrimaPagina)+"?s=t") # Se e' 3 va al trash
			
			
		zuppaPrimaPagina = BeautifulSoup(self.br.response().read())
	

		linkBase = self.funzLinkBase(zuppaPrimaPagina)
	
	# Cerca i messaggi nuovi: sono le table raw con sfondo bianco
	# e estrapola mittente e obj
		linksMailPrimaNew = zuppaPrimaPagina.findAll('tr', bgcolor="#ffffff")
		mittenteNew, oggettoNew, linkApriRelNew = self.estrapolaInfoNew(linksMailPrimaNew)
		linkApriNew = [linkBase+linkApriRelNew[i] for i in range(len(linkApriRelNew))]

	# Cerca i messaggi vecchi: sono le (table raw con) sfondo grigio
	# e estrapola mittente e obj
		linksMailPrimaOld = zuppaPrimaPagina.findAll(bgcolor="#E8EEF7")
		mittenteOld, oggettoOld, linkApriRelOld = self.estrapolaInfoOld(linksMailPrimaOld)
		linkApriOld = [linkBase+linkApriRelOld[i] for i in range(len(linkApriRelOld))]

		return mittenteNew, oggettoNew, linkApriNew, mittenteOld, oggettoOld, linkApriOld
# Fine funzione retriveMail(email, passwd)

# Funzione estrapolazione informazioni mail vecchie.
	def estrapolaInfoOld(self,links):
		mittente = []
		oggetto = []
		linkApriRel = []
		for i in range(0,len(links),1):
			str1 = str(links[i])
			zuppa = BeautifulSoup(str1)
			linkApriRelativo = zuppa.a.attrs[0][1]
			linkApriRel.append(linkApriRelativo)
			tdTrovati = zuppa.findAll('td')
		
			# Inizio estrapolazione mittente
			mittenteHtml = str(tdTrovati[1])
			mittenteTesto = mittenteHtml[5:len(tdTrovati[1])-6]
			mittente.append(mittenteTesto)
			# Fine estrapolazione mittente
		
			# Inizio estrapolazione oggetto
			spanTrovati = zuppa.findAll('span')
			oggettoHtml = str(spanTrovati[0])
			if oggettoHtml[55:68] =="Sent Messages":
				inizioOggettoTesto = 84
			else:
				inizioOggettoTesto = 70
				espressioneFine = 'font color="#7777CC"'
				espressioneFineComp = re.compile(espressioneFine)
				fineOggettoTestoS = re.search('font color="#7777CC"',oggettoHtml)

				if type(fineOggettoTestoS) == type(None):
					fineOggettoTesto = len(oggettoHtml)-8
				else:
					fineOggettoTesto = fineOggettoTestoS.start()-1	
			oggettoTesto = oggettoHtml[inizioOggettoTesto:fineOggettoTesto] 
			oggetto.append(oggettoTesto)
			# Fine estrapolazione oggetto
		
			# E' necessario ripulire il mittente (il primo messaggio "viene sporco")
		mittente = funzioneSpazzino.funzioneSpazzino(['idth="25%">', '</font>', '<font color="#777'], mittente)
		oggetto = funzioneSpazzino.funzioneSpazzino(['</font>', '<font color="#777', ", Sent Messages "], oggetto)
		return mittente, oggetto, linkApriRel
# Fine estrapolaInfoOld (links)	
	


# Funzione estrapolazione informazioni mail nuove.
	def estrapolaInfoNew(self,links):
		mittente = []
		oggetto = []
		linkApriRel = []
		for i in range(0,len(links),1):
			str1 = str(links[i])
			zuppa = BeautifulSoup(str1)
			linkApriRelativo = zuppa.a.attrs[0][1]
			linkApriRel.append(linkApriRelativo)
			tdTrovati = zuppa.findAll('td')

			# Inizio estrapolazione mittente
			mittenteHtml = str(tdTrovati[1])
			mittenteTesto = mittenteHtml[5:len(tdTrovati[1])-7]
			mittente.append(mittenteTesto)
			# Fine estrapolazione mittente
		
			# Inizio estrapolazione oggetto
			spanTrovati = zuppa.findAll('span')
			oggettoHtml = str(spanTrovati[0])
			if oggettoHtml[55:68] =="Sent Messages":
				inizioOggettoTesto = 87
			else:
				inizioOggettoTesto = 73
			espressioneFine = 'font color="#7777CC"'
			espressioneFineComp = re.compile(espressioneFine)
			fineOggettoTestoS = re.search('font color="#7777CC"',oggettoHtml)
			if type(fineOggettoTestoS) == type(None):
				fineOggettoTesto = len(oggettoHtml)-8
			else:
				fineOggettoTesto = fineOggettoTestoS.start()-6
			oggettoTesto = oggettoHtml[inizioOggettoTesto:fineOggettoTesto] 
			oggetto.append(oggettoTesto)
		# Fine estrapolazione oggetto
		
	# E' necessario ripulire il mittente
		mittente = funzioneSpazzino.funzioneSpazzino(['idth="25%">', '<b>','</b>','</t'], mittente)
		return mittente, oggetto, linkApriRel
	# Fine estrapolaInfoNew(links)

	def funzLinkBase(self,zuppaPrimaPagina):
		headPrima = zuppaPrimaPagina.findAll('head')
		linkBaseHtml = str(headPrima)
		linkBaseZuppa = BeautifulSoup(linkBaseHtml)
		return linkBaseZuppa.base.attrs[0][1]
		

# Prova apertura mail una per una
#	
#	def apriMail(self,link):
#		self.br.open(link)
#		zuppaMail = BeautifulSoup(self.br.response().read())
#		fileMail = open('zuppaMail.txt', 'w')
#		fileMail.write(str(zuppaMail))
#		fileMail.close()


