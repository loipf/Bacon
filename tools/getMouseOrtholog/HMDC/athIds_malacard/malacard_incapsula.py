#!/bin/python3
# incapsula cracker script to bypass crawler protection

import sys
from incapsula import IncapSession

session = IncapSession(user_agent='any-user-agent-string')
response = session.get('https://www.malacards.org/card/'+sys.argv[1])

print(response.text)








