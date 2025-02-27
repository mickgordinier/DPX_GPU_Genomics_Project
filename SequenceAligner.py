import numpy as np
from abc import ABC, abstractmethod

class SequenceAligner(ABC):
  
  def __init__(self, reference: str, query: str):
    self.reference = reference
    self.query = query
    self.Memo = None
  
  # Prints the matrix in a nice format
  def printMemoMatrix(self) -> str:
    lineToPrint = "      "
    for char in self.reference:
        lineToPrint += char + "  "
    print(lineToPrint)

    print(" ", self.Memo[0])
    for abc in range(1, len(self.query) + 1):
      print(self.query[abc - 1], self.Memo[abc])
  
  # Functions below must be overrided by the other sequence algorithms
  
  @abstractmethod
  def execute(self):
    pass
  
  @abstractmethod
  def initializeMemoMatrix(self):
    pass
  
  @abstractmethod
  def performRecursiveAnalysis(self):
    pass
  
  @abstractmethod
  def backtrack(self):
    pass
  