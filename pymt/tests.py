#!/usr/bin/env python

# Import libraries
import time
import hmac
import sha
import base64
import urllib
import xml.dom.minidom

# Define constants
AWS_ACCESS_KEY_ID = 'AKIAJJ66HQR34VAUE3DQ'
AWS_SECRET_ACCESS_KEY = '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'
SERVICE_NAME = 'AWSMechanicalTurkRequester'
SERVICE_VERSION = '2007-03-12'

sandbox = True
if sandbox:
  url = 'http://mechanicalturk.sandbox.amazonaws.com/onca/xml?'
else:
  url ='http://mechanicalturk.amazonaws.com/onca/xml?'

# Define authentication routines
def generate_timestamp(gmtime):
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", gmtime)

def generate_signature(service, operation, timestamp, secret_access_key):
    my_sha_hmac = hmac.new(secret_access_key, service + operation + timestamp, sha)
    my_b64_hmac_digest = base64.encodestring(my_sha_hmac.digest()).strip()
    return my_b64_hmac_digest

def acct():
  # Calculate the request authentication parameters
  operation = 'GetAccountBalance'
  timestamp = generate_timestamp(time.gmtime())
  signature = generate_signature('AWSMechanicalTurkRequester', operation, timestamp, AWS_SECRET_ACCESS_KEY)
  
  # Construct the request
  parameters = {
      'Service': SERVICE_NAME,
      'Version': SERVICE_VERSION,
      'AWSAccessKeyId': AWS_ACCESS_KEY_ID,
      'Timestamp': timestamp,
      'Signature': signature,
      'Operation': operation,
      }
  
  # Make the request
  result_xmlstr = urllib.urlopen(url, urllib.urlencode(parameters)).read()
  result_xml = xml.dom.minidom.parseString(result_xmlstr)
  
  # Check for and print results and errors
  errors_nodes = result_xml.getElementsByTagName('Errors')
  if errors_nodes:
      print 'There was an error processing your request:'
      for errors_node in errors_nodes:
  	for error_node in errors_node.getElementsByTagName('Error'):
  	    print '  Error code:    ' + error_node.getElementsByTagName('Code')[0].childNodes[0].data
  	    print '  Error message: ' + error_node.getElementsByTagName('Message')[0].childNodes[0].data
  
  availbalance_nodes = result_xml.getElementsByTagName('AvailableBalance')
  if availbalance_nodes:
      print "Available balance: " + availbalance_nodes[0].getElementsByTagName('FormattedPrice')[0].childNodes[0].data

def get_help():
  # Calculate the request authentication parameters
  operation = 'Help'
  timestamp = generate_timestamp(time.gmtime())
  signature = generate_signature('AWSMechanicalTurkRequester', operation, timestamp, AWS_SECRET_ACCESS_KEY)
  
  # Construct the request
  parameters = {
      'Service': SERVICE_NAME,
      'Version': SERVICE_VERSION,
      'AWSAccessKeyId': AWS_ACCESS_KEY_ID,
      'Timestamp': timestamp,
      'Signature': signature,
      'Operation': operation,
      }
  
  # Make the request
  result_xmlstr = urllib.urlopen(url, urllib.urlencode(parameters)).read()
  result_xml = xml.dom.minidom.parseString(result_xmlstr)
  
  # Check for and print results and errors
  errors_nodes = result_xml.getElementsByTagName('Errors')
  if errors_nodes:
      print 'There was an error processing your request:'
      for errors_node in errors_nodes:
  	for error_node in errors_node.getElementsByTagName('Error'):
  	    print '  Error code:    ' + error_node.getElementsByTagName('Code')[0].childNodes[0].data
  	    print '  Error message: ' + error_node.getElementsByTagName('Message')[0].childNodes[0].data
  
  print result_xml.toprettyxml()
  return result_xml

def getQuestion():
  return '''<QuestionForm xmlns="[the QuestionForm schema URL]">
  <Overview>
    <Title>Game 01523, "X" to play</Title>
    <Text>
      You are helping to decide the next move in a game of Tic-Tac-Toe.  The board looks like this:
    </Text>
    <Binary>
      <MimeType>
        <Type>image</Type>
        <SubType>gif</SubType>
      </MimeType>
      <DataURL>http://tictactoe.amazon.com/game/01523/board.gif</DataURL>
      <AltText>The game board, with "X" to move.</AltText>
    </Binary>
    <Text>
      Player "X" has the next move.
    </Text>
  </Overview>
  <Question>
    <QuestionIdentifier>nextmove</QuestionIdentifier>
    <DisplayName>The Next Move</DisplayName>
    <IsRequired>true</IsRequired>
    <QuestionContent>
      <Text>
        What are the coordinates of the best move for player "X" in this game?
      </Text>
    </QuestionContent>
    <AnswerSpecification>
      <FreeTextAnswer>
        <Constraints>
          <Length minLength="2" maxLength="2" />
        </Constraints>
        <DefaultText>C1</DefaultText>
      </FreeTextAnswer>
    </AnswerSpecification>
  </Question>
  <Question>
    <QuestionIdentifier>likelytowin</QuestionIdentifier>
    <DisplayName>The Next Move</DisplayName>
    <IsRequired>true</IsRequired>
    <QuestionContent>
      <Text>
        How likely is it that player "X" will win this game?
      </Text>
    </QuestionContent>
    <AnswerSpecification>
      <SelectionAnswer>
        <StyleSuggestion>radiobutton</StyleSuggestion>
        <Selections>
          <Selection>
            <SelectionIdentifier>notlikely</SelectionIdentifier>
            <Text>Not likely</Text>
          </Selection>
          <Selection>
            <SelectionIdentifier>unsure</SelectionIdentifier>
            <Text>It could go either way</Text>
          </Selection>
          <Selection>
            <SelectionIdentifier>likely</SelectionIdentifier>
            <Text>Likely</Text>
          </Selection>
        </Selections>
      </SelectionAnswer>
    </AnswerSpecification>
  </Question>
</QuestionForm>'''

def make2():
  import boto
  mt = boto.connect_mturk( AWS_ACCESS_KEY_ID,  AWS_SECRET_ACCESS_KEY)
  return mt

def makeHit():
  # Calculate the request authentication parameters
  operation='CreateHIT'
  timestamp = generate_timestamp(time.gmtime())
  signature = generate_signature('AWSMechanicalTurkRequester', operation, timestamp, AWS_SECRET_ACCESS_KEY)
  

  HITTypeId='T100CN9P324W00EXAMPLE'
  question=[getQuestion()]
  lifetimeInSeconds=604800
  assignmentDurationInSeconds=100
  reward = '''<Reward>
  <Amount>0.05</Amount>
  <CurrencyCode>USD</CurrencyCode>
</Reward>'''
          

  # Construct the request
  parameters = {
      'Service': SERVICE_NAME,
      'Version': SERVICE_VERSION,
      'AWSAccessKeyId': AWS_ACCESS_KEY_ID,
      'Timestamp': timestamp,
      'Signature': signature,
      'Operation': operation,
      'LifetimeInSeconds': lifetimeInSeconds,
      'Question': question,
      'Title':'This is a test hit',
      'Description':'This hit is a demo from the website. Do what you like with it.',
      'AssignmentDurationInSeconds':assignmentDurationInSeconds,
      'Reward':reward
      }
  
  # Make the request
  result_xmlstr = urllib.urlopen(url, urllib.urlencode(parameters)).read()
  result_xml = xml.dom.minidom.parseString(result_xmlstr)
  
  # Check for and print results and errors
  errors_nodes = result_xml.getElementsByTagName('Errors')
  if errors_nodes:
      print 'There was an error processing your request:'
      for errors_node in errors_nodes:
  	for error_node in errors_node.getElementsByTagName('Error'):
  	    print '  Error code:    ' + error_node.getElementsByTagName('Code')[0].childNodes[0].data
  	    print '  Error message: ' + error_node.getElementsByTagName('Message')[0].childNodes[0].data
  
  print result_xml.toprettyxml()
  return result_xml

if __name__ == '__main__':
  import sys
  args = sys.argv[1:]
  if len(args) == 0:
    get_help()
  else:
    if args[0] == 'balance': acct()
    print 'sorry no args are yet implemented'
