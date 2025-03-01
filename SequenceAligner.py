import numpy as np
from abc import ABC, abstractmethod

RED = '\033[31m'
GREEN = '\033[32m'
RESET = '\033[0m'

class SequenceAligner(ABC):
  
  def __init__(self, reference: str, query: str):
    self.reference = reference
    self.query = query
    self.Memo = None
  
  # Prints the matrix in a nice format
  def printMemoMatrix(self) -> str:
    lineToPrint = "          "
    for char in self.reference:
        lineToPrint += char + "    "
    print(RED + lineToPrint + RESET)
      
    # Print numpy matrix using the numpy matrix print and save it to a string
    # Pad each digit to 3 decimal places for nice formatting
    matrix_str = np.array2string(self.Memo, 
                             formatter={'all': lambda x: f'{int(x):3d}'}, 
                             separator=', ', 
                             max_line_width=80, 
                             precision=0)
    
    # Remove the beginning and end backets
    matrix_str = matrix_str[1:-1]
    
    # Remove the spacing after each newline representing a new row
    matrix_rows = [row.strip() for row in matrix_str.split("\n")]
    
    # Add each character of query to each row
    matrix_rows[0] = "  " + matrix_rows[0]
    for row_idx in range(1, len(self.query)+1):
      matrix_rows[row_idx] = f"{GREEN}{self.reference[row_idx]}{RESET} {matrix_rows[row_idx]}"
    
    matrix_str = "\n".join(matrix_rows)
    
    print(matrix_str)
  
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
  