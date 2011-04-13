* * * * * *

Getting Started Sample for Amazon Mechanical Turk, in Python (REST)

* * * * * *

An example of making a call to the Amazon Mechanical Turk web service,
using the REST interface, in Python.


About This Sample

 * This is an example of making a call to the Amazon Mechanical Turk
web service, using the REST interface, in Python.

 * This example demonstrates how to sign an AWS request in Python. 

 * This code is pre-built to use API version 2007-03-12.


Prerequisites

 * You will need an Amazon Mechanical Turk Requester account.  You can
sign up at the Requester web site: http://requester.mturk.amazon.com

 * You will to register for Amazon Web Services.  Be sure to use the
same e-mail address and password as you used when creating your
Requester account.  You can sign up at the AWS web site:
http://aws.amazon.com/register

 * This Python example uses Python 2.4.


Using the Sample

 1. Update the DoGetAccountBalance.py file to include your Amazon Web
Services Access Key ID and Secret Access Key.  You can retrieve this
information from the AWS web site: http://aws.amazon.com

 2. Run the example by entering the following command:

   python DoGetAccountBalance.py

The command response should be similar to the following:

   Account balance: $0.00
