"""### passwordGetter ###

usage: user, pass = passwordGetter()

The function passwordGetter() displays a simple Tkinter interface where the first entry
gets the User and the second gets the password.  
Clicking the button labeled 'Log In' the interface will close and the function will return
the values of the user-field and of the password-field.


Claudio Virili     claudio.virili@gmail.com
Emanuele Casali    emanuele.casali@gmail.com

22/01/2010
"""

from Tkinter import *

class MyApp:
    def __init__(self, parent):
	self.root1 = parent
        self.myContainer1 = Frame(parent)
        self.myContainer1.pack()
	Label(self.myContainer1, text = "User").grid(row=0,column=0)
        self.user_entry = Entry(self.myContainer1)
        self.user_entry.grid(row=0,column=1)

		

	Label(self.myContainer1, text = "Password").grid(row=1,column=0)
        self.password_entry = Entry(self.myContainer1, show="*")
        self.password_entry.grid(row=1,column=1)	 


	Button(self.myContainer1, text="Log In", width=10, command=self.callback).grid(row=2, column=1)
	  

    def callback(self):
	self.password = self.password_entry.get()   
	self.user = self.user_entry.get()
	self.root1.destroy()

def passwordGetter():
	root = Tk()
	myapp = MyApp(root)
	root.mainloop()
	return myapp.user, myapp.password
	
