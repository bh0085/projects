"""### gmail interface (main) ### 

Claudio Virili     claudio.virili@gmail.com
Emanuele Casali    emanuele.casali@gmail.com

22/01/2010
"""

from Tkinter import *
import gmailIM
import passwordnew


class finestraMsg:
	def __init__(self, email, password, parent):
		self.connessioneGmail = gmailIM.retriveMail(email, password)
		
		self.pagina = 1
		self.maxPagina = 5 # Per ora e' fisso
		
		self.MainFrame = Frame(parent, width=77, height=22)
		self.MainFrame.pack(side=TOP)
		
		# UPPER FRAME
		self.UpperFrame = Frame(self.MainFrame)
		self.parametroPagina = "1"
		self.frameAdattivo()
		
		# LOWER FRAME
		LowerFrame = Frame(parent)
		Button(LowerFrame, text="<<", command=self.decrementaPagina).grid(row=0,column=0)
		Button(LowerFrame, text="Inbox", command=lambda parametroPagina = "1": self.chiamataBottone(parametroPagina)).grid(row=0,column=1)
		Button(LowerFrame, text="Sent", command=lambda parametroPagina = "2": self.chiamataBottone(parametroPagina)).grid(row=0,column=2)
		Button(LowerFrame, text="Trash", command=lambda parametroPagina = "3": self.chiamataBottone(parametroPagina)).grid(row=0,column=3)
		Button(LowerFrame, text=">>", command=self.incrementaPagina).grid(row=0,column=4)
		LowerFrame.pack(side=BOTTOM)
	
	def incrementaPagina(self):
		if self.pagina != 5:
			self.pagina = self.pagina+1
		print self.pagina
		self.frameAdattivo()
			
	def decrementaPagina(self):
		if self.pagina != 1:
			self.pagina = self.pagina-1
		print self.pagina
		self.frameAdattivo()
		
	def chiamataBottone(self,parametroCiccio):
		self.pagina = 1
		self.parametroPagina=parametroCiccio
		self.frameAdattivo()
	
	def frameAdattivo(self):
		
		mittenteNew, oggettoNew, linkApriNew, mittenteOld, oggettoOld, linkApriOld = self.connessioneGmail.retriveMainInfos(self.parametroPagina)
		
		self.UpperFrame.destroy() 
		self.UpperFrame = Frame(self.MainFrame)
		
		# UPPER FRAME
		LabelTextmittente = mittenteNew+mittenteOld
		LabelTextoggetto = oggettoNew+oggettoOld
		UpperLabelmittente = []
		UpperLabeloggetto = []
		#BottoneAzione = []
	
		# Prima riga del frame = intestazioni
		UpperLabelmittente.append(Label(self.UpperFrame, text="Mittente", bg='red', width=25, anchor=W, relief=RIDGE).grid(row=0,column=0))
		UpperLabeloggetto.append(Label(self.UpperFrame, text="Oggetto", bg='red', width=50, anchor=W, relief=RIDGE).grid(row=0,column=1))
		#BottoneAzione.append(Label(UpperFrame, text="Azione", bg='red', width=10, anchor=W, relief=RIDGE).grid(row=0,column=2))
		# Ad ogni riga viene assegnata una mail 
		# Le mail nuove hanno sfondo bianco, le mail vecchie hanno sfondo grigio (come su gmaill.com)
		for i in range(self.pagina*10-10,self.pagina*10,1):
			if i < len(mittenteNew):
				UpperLabelmittente.append(Label(self.UpperFrame, text=LabelTextmittente[i], bg='white', width=25, height=2, anchor=W, relief=RIDGE).grid(row=i+1,column=0))
				UpperLabeloggetto.append(Label(self.UpperFrame, text=LabelTextoggetto[i], bg='white', width=50,  height=2, anchor=W, relief=RIDGE).grid(row=i+1,column=1))	
			else:
				UpperLabelmittente.append(Label(self.UpperFrame, text=LabelTextmittente[i], bg='grey', width=25, height=2, anchor=W, relief=RIDGE).grid(row=i+1,column=0))
				UpperLabeloggetto.append(Label(self.UpperFrame, text=LabelTextoggetto[i], bg='grey', width=50, height=2, anchor=W, relief=RIDGE).grid(row=i+1,column=1))	
				#BottoneAzione.append(Button(UpperFrame, text="Apri", width=5, height=1, anchor=W).grid(row=i+1,column=2))
	
		self.UpperFrame.pack()
	

	
		
email, password = passwordnew.passwordGetter()
MainWindow = Tk()
finestraMsg(email, password, MainWindow)
MainWindow.mainloop()