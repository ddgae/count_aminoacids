#Example of binding a sequence tag with protein/nucleic acid sequence. 
#I will continue to create sample code as tutorial for manupulating sequences. 
import os,sys,zip 

#example2 (simple solution)
a1 = list()
a2 = list()
a3 = []
a4 = []
total = list()
for x in seq:
  a1.append(x)
for y in quality:
  a2.append(y)

for i, (a3, a4) in enumerate(zip(a1,a2)):
  total = [name[i],[a3,a4]]
  print(total)


#example3 ( linked list solution)
class ListNode:
   #initial node
   def __init__(self, data=None, next=None):
       self.data = data
       self.next = next
   #return the node
   def __init__(self):
       return init(self.data)

class LinkedList:
   #start of the node
   def __init__(self):
       self.head = None
   #start at the node and add.
   def attend(self, data):
       self.head = ListNode(data=data, next=self.head)
if __name__ == "__main__":
   ll = LinkedList()
   ll.attend('seq_1    gtacgacgatcgactagc wiuehrwiuyeiuwh')
   ll.attend("seq_2    gacgtacgtacgactacgtcga kjfy3wrfhiwgiyfwd")
   ll.attend("seq_3    nnnnnnnnnnnnnn  jwerhfiuwijhwrfksjd")
