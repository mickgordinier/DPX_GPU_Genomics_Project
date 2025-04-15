from playwright.sync_api import sync_playwright
import time


def convert_to_dna(sequence):
    # Mapping from "0123" to "ATCG"
    mapping = {'0': 'A', '1': 'T', '2': 'C', '3': 'G', '4': 'U', '_': '_'}
    
    # Convert the sequence using the mapping
    converted_sequence = ''.join(mapping[digit] for digit in sequence)
    
    return converted_sequence


def readInputSequences(file_path):
  
  input_references = []
  input_queries = []
  
  with open(file_path, 'r') as file:
    for index, line in enumerate(file):
      
      if (index % 3 == 1): 
        # print("NEW REFERENCE:", line.strip())
        input_references.append(convert_to_dna(line.strip()))
      if (index % 3 == 2): 
        # print("NEW QUERY:", line.strip())
        input_queries.append(convert_to_dna(line.strip()))
      
  return input_references, input_queries


def readOutputSequences(filePath):
  scoreOutputs = []
  referenceOutputs = []
  alignmentOutputs = []
  queryOutputs = []

  with open(filePath, "r") as f:
      lines = f.readlines()

  lines = lines[0:-1]

  # Process every 4 lines in the remaining list
  for i in range(0, len(lines), 4):
    scoreOutputs.append(lines[i].split("|")[1].strip())
    referenceOutputs.append(convert_to_dna(lines[i + 1].strip()))
    alignmentOutputs.append(lines[i + 2].strip())
    queryOutputs.append(convert_to_dna(lines[i + 3].strip()))
    
  return scoreOutputs, referenceOutputs, alignmentOutputs, queryOutputs


def scrape_needleman_wunsch(page, reference, query):
    
  page.wait_for_selector("#sequence_1")
  page.wait_for_selector("#sequence_2")

  # Fill in the sequences
  page.fill("#sequence_1", reference)
  page.fill("#sequence_2", query)
  
  page.wait_for_selector("#match")
  page.fill("#match", "3")
  
  time.sleep(5)

  # Wait for the results table to be present
  page.wait_for_selector("table.results")

  # Find all rows in the results table
  rows = page.locator("table.results tr")
  
  results_output = []
  
  # Wait for at least one alignment to be present before continuing
  page.wait_for_selector('[data-bind^="text: $root.output.alignments"]', timeout=10000)

  # Iterate through each row and get the data bound to the text
  for i in range(rows.count()):
      row = rows.nth(i)
      
      # Use try-except to catch timeouts and handle them gracefully
      try:
        # Wait for each alignment text to exist inside the row before trying to read it
        cell1 = row.locator('[data-bind="text: $root.output.alignments()[$index()][0]"]')
        cell2 = row.locator('[data-bind="text: $root.output.alignments()[$index()][1]"]')
        cell3 = row.locator('[data-bind="text: $root.output.alignments()[$index()][2]"]')

        cell1.wait_for(timeout=5000)
        cell2.wait_for(timeout=5000)
        cell3.wait_for(timeout=5000)

        text1 = cell1.text_content(timeout=5000) or ""
        text2 = cell2.text_content(timeout=5000) or ""
        text3 = cell3.text_content(timeout=5000) or ""

        results_output.append([text1.strip(), text2.strip(), text3.strip()])
      
      except Exception as e:
          print(f"Row {i+1} - Error retrieving data: {e}")
          
  page.wait_for_selector('[data-bind="text: $root.output.score"]', timeout=10000)  # Wait up to 10 seconds
  
  # Locate the element with the specific data-bind
  score_element = page.locator('[data-bind="text: $root.output.score"]')
  
  # Extract the text content from the located element
  score_text = score_element.text_content(timeout=5000)  # Timeout for getting text
  
  return score_text, results_output
        
        
def main():
  
  print("READING INPUT FILES\n")
  
  input_references, input_queries = readInputSequences("input-data.txt")
  scoreOutputs, referenceOutputs, alignmentOutputs, queryOutputs = readOutputSequences("align-output.txt")
  
  print("FINISHED READING INPUT FILES\n")
  
  with sync_playwright() as p:
    
    # Launch Chromium browser
    browser = p.chromium.launch(headless=False)  # Set headless=True for invisible browser
    page = browser.new_page()

    # Open the target URL
    page.goto("https://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Gotoh")

    # Optional: Wait for page to load fully (e.g., wait for a specific element)
    page.wait_for_load_state("networkidle")
    
    page.wait_for_selector("#match")
    page.wait_for_selector("#mismatch")
    page.wait_for_selector("#base_costs")
    page.wait_for_selector("#enlargement")
    
    page.fill("#match", "3")
    page.fill("#mismatch", "-1")
    page.fill("#base_costs", "-3")
    page.fill("#enlargement", "-1")
    
    for i in range(len(input_references)):
    
      print("------------------------------- PERFORMING NEW ALIGNMENT TESTING -----------------------------------\n")
      print("PAIR #:", i)
      print("INPUT REFERENCE:", input_references[i])
      print("INPUT QUERY:", input_queries[i])
      print()
      
      score, paths = scrape_needleman_wunsch(page, input_references[i], input_queries[i])
      
      print("*** CALCULATED SCORES ***")
      print("SCORE:", scoreOutputs[i])
      print("REFERENCE:", referenceOutputs[i])
      print("ALIGNMENT:", alignmentOutputs[i])
      print("QUERY:    ", queryOutputs[i])
      print()
      
      print("*** OUTPUT SCORES ***")
      print("SCORE:", score)
      
      if (score != scoreOutputs[i]):
        print("ERROR: SCORES DON'T MATCH!!!")
        exit(1)
      
      print("SCORES MATCH!!!\n")
      
      if (len(paths) == 10):
        print("There exists too many path to compare")
        print("----------------------------------------------------------------------------------------------------\n\n")
        continue
      
      print("There are less than 10 paths, can view them all!\n")
      
      pathFound = False
      
      for onlinePath in paths:
        if ((onlinePath[0] == referenceOutputs[i]) and (onlinePath[1] == alignmentOutputs[i]) and (onlinePath[2] == queryOutputs[i])):
          print("FOUND!!!!")
          pathFound = True
        else:
          print("Different online values")
        
        print("REFERENCE:", onlinePath[0])
        print("ALIGNMENT:", onlinePath[1])
        print("QUERY:    ", onlinePath[2])
        print()
          
      if (not pathFound):
        print("ERROR: MATCHING PATH COULD NOT BE FOUND!!!")
        exit(1)
          
      print("----------------------------------------------------------------------------------------------------\n\n")
    
      
    # Close browser
    browser.close()
    
  #end with
  
  

if __name__ == "__main__":
    main()