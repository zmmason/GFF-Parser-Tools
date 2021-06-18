# !/usr/bin/env python3
# University of California, Santa Cruz
# Biomolecular Engineering and Bioinformatics
# Name: Zachary Mason (zmmason), Justin Chan (jumchan)
# Group Members: None


class ProteinParam:
    """
    The class, ProteinParam, is formatted to have one constructor and seven
    methods that calculate their own quantitative characteristics of the 
    input amino acid sequence.
    """
    aa2mw = {
        'A': 89.093, 'G': 75.067, 'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
    }

    mwH2O = 18.015
    aa2abs280 = {'Y': 1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R': 12.4, 'H': 6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__(self, protein):
        """
        To make a ProteinParam object, a input string that represents an amino
        acid sequence need to be passed in.
        
        The input sequence is then formatted in all uppercase and valid single
        letter amino acid codes in a string and then formatted with a dictionary
        of values with the counts of every valid amino acid.
        
        This design of a dictionary with the count of amino acids was decided
        because it would be the easiest to have when calculating return values
        for the methods, molarExtinction(), molecularWeight(), and charge().
        """
        realAASeq = ""
        listOfChar = list(protein.upper())  # Makes AA sequence input capitalized
        realAA = self.aa2mw.keys()  # Creates list of keys for valid amino acids
        for aminoAcid in listOfChar:  # Checks to see if all characters are valid amino acid codes
            for realAminoAcid in realAA:
                if aminoAcid == realAminoAcid:
                    realAASeq += aminoAcid

        aaSeq = realAASeq

        aminoAA = list(self.aa2mw.keys())
        numOfAA = [0] * len(aminoAA)  # Creates beginning counts of each valid amino acid
        indivAAList = list(aaSeq)
        for aminoAcid in indivAAList:  # Counts how many of each valid amino acids there are
            for singleAA in aminoAA:
                if singleAA == aminoAcid:
                    numOfAA[aminoAA.index(singleAA)] += 1

        self.countAA = dict(zip(aminoAA, numOfAA))  # Creates a dictionary with amino acid counts
        # Single letter amino acid codes are keys and their counts in the input sequence are their values

    def aa_Count(self):
        """
        Returns how total valid single letter amino acid codes in the input amino acid sequence 
        
        Calculated by using a for loop through the dictionary of amino acid counts to return the
        sum of all the counts
        
        input: dictionary containing amino acid counts of input sequence
        output: valid number of amino acids in input sequence
        """
        count = 0
        for countAA in self.countAA.values():  # Counts total amino acids by taking sum of counts in dictionary
            count += countAA
        return count

    def pI(self):
        """
        Returns the theoretical isoelectric point. The theoretical isoelectric point is the pH
        where the input amino acid sequence is neutralized or when its net charge is closest to 0
        
        The method pI is basically a for loop that loops through pHs 0-14 by every 0.01 pH and
        uses the charge() to check the charge at that pH. If the charge is lower than the saved
        charge, it saves the charge and that pH. This is done to decide which pH provides the
        most neutral charge for the amino acid sequence and returns the pH or theoretical 
        isoelectric point.
        
        input: dictionary containing amino acid counts of input sequence
        output: the pH where the amino acid sequence is neutralized or has a net charge closest to 0
        """
        bestCharge = 1000000  # Represents smallest charge through the pH range 0-14
        bestpH = 0.00  # Represents pH with the most neutral chrage for the amino acid sequence
        for pH in range(0, 1401, 1):  # Loops through pH 0-14 every 0.01
            realpH = pH / 100
            thisCharge = abs(self.charge(realpH))  # Absolute value of charge at pH using charge()
            if thisCharge < bestCharge:  # Checks if the recent pH produces a smaller charge
                bestCharge = thisCharge
                bestpH = realpH
        return bestpH

    def aaComposition(self):
        """
        Return a dictionary containing valid single letter amino acid codes and their counts. The method
        is pretty simple because the dictionary is made in the constructor and aaComposition() simply
        returns it.
        
        input: dictionary containing amino acid counts of input sequence
        output: dictionary containing amino acid counts of input sequence
        """
        return self.countAA

    def charge(self, pH):
        """
        Returns the charge of the amino acid sequence given a pH. 
        
        In other words, charge returns the charge of an amino acid sequence given a passed in pH.
        The method charge() is private because its sole use is only to be used in another method, pI().
        
        charge() utilizes a for loop and a series of if statements to determine if amino acid
        is negatively and positively charged in the amino acid dictionary of counts given the dictionaries
        aa2chargePos and aa2chargeNeg. If it is, it uses an equation with the amino acid's pKa and current 
        pH that adds to the positive or negative charge to be used to calculate the net charge. It also 
        factors the N and C terminus into the charge() after using object attributes, aaNterm and aaCterm
        and return the net charge.
        
        Positvely Charged AAs: K,R,H (N terminus)
        Negatively Charged AAs: D,E,C,Y (C terminus)
        
        Class Attributes Used: aa2chargePos, aa2chargeNeg, aaNterm, aaCterm
        The class attributes used were need to get the pKa's of the amino acids which were saved in 
        floats and in dictionaries in the class to find the net charge of the protein sequence.
        
        input: dictionary containing amino acid counts of input sequence
        output: charge of the amino acid at a specific pH
        """
        posCharge = 0.00  # Charge contribution from positively charged amino acids and N terminus
        negCharge = 0.00  # Charge contribution from negatively charged amino acids and C terminus
        if self.aa_Count() != 0:
            for aminoAcid, countAminoAcid in self.countAA.items():  # Loops through amino acids
                if self.aa2chargePos.get(aminoAcid,
                                         'unknown') != 'unknown':  # Checks for positively charged amino acids
                    posCharge += countAminoAcid * (10 ** (self.aa2chargePos.get(aminoAcid))) / (
                            10 ** self.aa2chargePos.get(aminoAcid) + 10 ** pH)
                elif self.aa2chargeNeg.get(aminoAcid,
                                           'unknown') != 'unknown':  # Checks for negatively charged amino acids
                    negCharge += countAminoAcid * (10 ** pH) / (10 ** self.aa2chargeNeg.get(aminoAcid) + 10 ** pH)

            posCharge += (10 ** self.aaNterm) / (10 ** self.aaNterm + 10 ** pH)  # Factors N terminus in net charge
            negCharge += (10 ** pH) / (10 ** self.aaCterm + 10 ** pH)  # Factors C terminus in net charge

        return posCharge - negCharge

    def molarExtinction(self):
        """
        Returns the molar extinction coefficent. The molar extinction coefficient is how much light
        a protein absorbs at a given wavelength.
        
        A for loop is used to go through every item of the amino acid count dictionary to determine
        if an amino acid is one that is part of the calculation of the molar extinction coefficient.
        If the amino acid is either Y, W, C, then it is factored in calculating the molar extinction
        coefficient by multiplying its count and extinct coefficients.
        
        Class Attributes Used: aa2abs280
        The class attribute used here was used to find the molar extinction coefficients of particular
        amino acids that contributed to the molar extinction coefficient of the total protein sequence.
        
        input: dictionary containing amino acid counts of input sequence
        output: molar extinction coefficient of a protein
        """
        molEx = 0.000  # Molar Extinction Coefficient of AA Sequence
        for aminoAcidAA, countAminoAcidAA in self.countAA.items():  # Loops through AA count dictionary
            if self.aa2abs280.get(aminoAcidAA, 'unknown') != 'unknown':  # Checks if amino acid affects light absorption
                molEx += countAminoAcidAA * self.aa2abs280.get(aminoAcidAA)
        return molEx

    def massExtinction(self):
        """
        Returns mass extinction coefficient. Divides molar extinction coefficient by molecular weight
        to calculate mass extinction coefficient
        
        massExtinction() utilizes molecularWeight() and molarExtinction() methods to return mass 
        extinction coefficient. 
        
        input: dictionary containing amino acid counts of input sequence
        output: mass extinction coefficient of a protein
        """
        myMW = self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0  # Ensures molecular weight is greater than 0

    def molecularWeight(self):
        """
        Returns the molecular weight of the amino acid sequence by adding weight of amino acid
        sequence.
        
        molecularweight() goes through the amino acid count dictionary and counts how many of each
        amino acid there is and uses it to multiply it by its molecular weight in the aa2mw dictionary.
        It simply then returns sum of these amino acids in the sequence.
        
        Class Attributes Used: aa2mw
        This attribute was used to determine the molecular weight of each amino acid saved as
        a dictionary.
        
        input: dictionary containing amino acid counts of input sequence
        output: molecular weight of the amino acid sequence
        """
        weight = 0.000
        if self.aa_Count() != 0:
            for aminoAA, countAA in self.countAA.items():  # Loops through each dictionary item of countAA
                weight += self.aa2mw.get(aminoAA) * countAA
                weight -= self.mwH2O * countAA  # Get rid of H20 weight in the amino acid
            weight += self.mwH2O
        return weight


class NucParams:
    """
    The class, NucParams, has one constructor and five methods that use to count amino acid,
    single nucleotides, and codons of an DNA/RNA nucleotide sequence.
    
    The class has 3 class attributes that are dictionaries that relate counts of nucleotides
    and codons in the input string and the relationships between codons and amino acid letter
    codes. These dictionaries are used in all the methods to return counts and the dictionaries
    themselves to be used to gain information about a DNA/RNA nucleotide sequence.
    
    NucParams takes an input string and runs for loops that stores single nucleotides and codons.
    This was done because given the methods, it would be much easier to create dictionaries of 
    the input sequence in the constructor using addSequence() and return them in methods. Furthermore,
    many of the methods rely a lot on each other so it seemed best to make them attributes of the class.
    
    Three of the methods in the class return dictionaries of valid nucleotides, codons, and amino acids
    in the input sequnece. One method returns the total count of nucleotides and the last method allows
    to add on the previous sequence with new sequence that can be analyzed together.
    """

    rnaCodonTable = {
        # U
        'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
        'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
        'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
        'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
        # C
        'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
        'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
        'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
        'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
        # A
        'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
        'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
        'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
        'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
        # G
        'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
        'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
        'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
        'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U', 'T'): value for key, value in rnaCodonTable.items()}

    def __init__(self, inString=''):
        """
        The constructor takes in an optional string parameter that is suppose to represent the DNA/RNA
        nucleotide sequence. If there is no sequence, it gives a default empty string.
        
        The constructor saves two additional dictionaries to the class that count single nucleotides 
        and codons of the input string.
        
        To get these counts, the constructor uses addSequence() to first capitalize all letters and 
        gets rid of spaces before adding to dictionaries. Then it uses for loops to go through the 
        input string by character to determine if it is a valid single nucleotide and then by a string
        of three characters to determine of it a valid codon and saves these in dictionaries.
        
        The design of dictionaries was used because it seemed that since the constructor and addSequence()
        do practically the same thing, sort nucleotides and codons into dictionaries and many of these
        methods involve returning these dictionaries, it was decided that it was better to have class
        attributes as dictionaries.
        
        """
        self.validNucDic = {'A': 0, 'C': 0, 'T': 0, 'G': 0, 'U': 0, 'N': 0}  # Valid nucleotide count dictionary
        self.rnaCodonCount = {
            # Valid codon count dictionary
            # U
            'UUU': 0, 'UCU': 0, 'UAU': 0, 'UGU': 0,  # UxU
            'UUC': 0, 'UCC': 0, 'UAC': 0, 'UGC': 0,  # UxC
            'UUA': 0, 'UCA': 0, 'UAA': 0, 'UGA': 0,  # UxA
            'UUG': 0, 'UCG': 0, 'UAG': 0, 'UGG': 0,  # UxG
            # C
            'CUU': 0, 'CCU': 0, 'CAU': 0, 'CGU': 0,  # CxU
            'CUC': 0, 'CCC': 0, 'CAC': 0, 'CGC': 0,  # CxC
            'CUA': 0, 'CCA': 0, 'CAA': 0, 'CGA': 0,  # CxA
            'CUG': 0, 'CCG': 0, 'CAG': 0, 'CGG': 0,  # CxG
            # A
            'AUU': 0, 'ACU': 0, 'AAU': 0, 'AGU': 0,  # AxU
            'AUC': 0, 'ACC': 0, 'AAC': 0, 'AGC': 0,  # AxC
            'AUA': 0, 'ACA': 0, 'AAA': 0, 'AGA': 0,  # AxA
            'AUG': 0, 'ACG': 0, 'AAG': 0, 'AGG': 0,  # AxG
            # G
            'GUU': 0, 'GCU': 0, 'GAU': 0, 'GGU': 0,  # GxU
            'GUC': 0, 'GCC': 0, 'GAC': 0, 'GGC': 0,  # GxC
            'GUA': 0, 'GCA': 0, 'GAA': 0, 'GGA': 0,  # GxA
            'GUG': 0, 'GCG': 0, 'GAG': 0, 'GGG': 0  # GxG
        }

        self.addSequence(inString)  # Calls addSequence() to add to the dictionary counts

    def addSequence(self, inSeq):
        """
        Adds a new DNA/RNA nucleotide sequence to the previous DNA/RNA nucleotide sequence
        to be counted.
        
        An input string is taken as parameter which represents the sequence added to the previous
        sequence.
        
        The method is used to format the string and then use it to add to the counts of valid
        nucleotides and codons. addSequence() capaitalizes and gets rid of white spaces to format
        the DNA/RNA sequence. It then loops through the string by charcater and then by three letters
        to check for valid nucleotides and codons and adds them to the dictionary.
        
        Class Attributes Used: validNucDic, rnaCodonCount
        
        input: an input string that represents a DNA/RNA sequence
        output: no output 
        """
        addedDNACap = inSeq.upper()  # Capitalizes DNA/RNA sequence
        addedDNANoSpace = addedDNACap.replace(' ', '')  # Gets rid of spaces
        for singNuc in addedDNANoSpace:
            if singNuc in self.validNucDic:
                self.validNucDic[singNuc] += 1  # Adds to count
        rnaSeqConv = addedDNANoSpace.replace('T', 'U')
        for i in range(0, len(rnaSeqConv), 3):
            if rnaSeqConv[i:i + 3] in self.rnaCodonCount:
                self.rnaCodonCount[rnaSeqConv[i:i + 3]] += 1  # adds to count

    def aaComposition(self):
        """
        Determines the amino acid counts of the 20 amino acids and stop codon and
        returns the counts in a dictionary.
        
        Creates a dictionary called aa_Count which is an dictionary with all amino acids
        that begin with a count 0. The codon count dictionary is then looped and checks
        to see which codon belongs to which amino acid using rnaCodonTable and will add 
        to a item in aa_Count if it determines a codon belongs to a particular amino acid.
        
        Class Attributed Used: rnaCodonCount, rnaCodonTable
        
        input: A RNA codon count dictionary that counts codons in a DNA/RNA sequence
        output: amino acid count dictionary in a RNA/DNA sequence
        """
        aa_Count = {
            'F': 0, 'L': 0, 'I': 0, 'M': 0, 'V': 0,
            'S': 0, 'P': 0, 'T': 0, 'A': 0, 'Y': 0,
            'H': 0, 'Q': 0, 'N': 0, 'K': 0, 'D': 0,
            'E': 0, 'C': 0, 'W': 0, 'R': 0, 'G': 0, '-': 0
        }
        for codonCode in self.rnaCodonCount.keys():
            if self.rnaCodonCount.get(codonCode) != 0:  # Checks if count is not 0
                aa_Count[self.rnaCodonTable.get(codonCode)] += self.rnaCodonCount.get(codonCode)
        return aa_Count

    def nucComposition(self):
        """
        Returns a dictionary containing valid single nucleotide counts in the DNA/RNA
        sequence.
        
        The method is very simple which it simply returns the class attribute of
        validNucDic
        
        Class Attributes Used: validNucDic
        
        input: A dictionary containing single nucleotide counts in input sequence
        output: A dictionary containing single nucleotide counts in input sequence
        """
        return self.validNucDic

    def codonComposition(self):
        """
        Returns a dictionary containing valid rna codon counts in the DNA/RNA
        sequence.
        
        The method is very simple which it simply returns the class attribute of
        rnaCodonCount
        
        Class Attributes Used: rnaCodonCount
        
        input: A dictionary containing valid codon counts in input sequence
        output: A dictionary containing valid codon counts in input sequence
        """
        return self.rnaCodonCount

    def nucCount(self):
        """
        Returns a count of total valid nucleotides in the DNA/RNA sequence.
        
        The method is very simple which it simply returns an integer of total
        nucleotides in the sequence.
        
        Class Attributes Used: validNucDic
        
        input: A dictionary containing valid single nucleotide counts in input sequence.
        output: An integer of how many total valid nucleotides there are in the sequence.
        """
        count = 0
        for nucC in self.validNucDic.values():  # Loops through nucleotide counts in validNucDict
            count += nucC  # adds counts to total count
        return count


import sys


class FastAreader:
    """ 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    """

    def __init__(self, fname=''):
        """contructor: saves attribute fname """
        self.fname = fname

    def doOpen(self):
        """ Handle file opens, allowing STDIN."""
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):
        """ Read an entire FastA record and return the sequence header/sequence"""
        header = ''
        sequence = ''
        with self.doOpen() as fileH:
            header = ''
            sequence = ''
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
            header = line[1:].rstrip()
            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()
        yield header, sequence


class OrfFinder():
    """
    Create a list of open reading frames that contain its frame, positions, and length
    
    The OrfFinder is set up to add open reading frames to openReadingFrames, a class 
    attribute that is a list of all the ORFs found in a DNA sequence. openReadingFrames
    is needed to store the ORFs found in the sequence and is used to create the output
    file.
    
    The class contains a constructor and three methods that are used to generate ORFs
    and to add them to openReadingFrames. A majority of the program is done in findOrfs() 
    because it felt it would be simplier to find and keep positons, frames, and length of ORFs
    done in one method which could be easily added to openReadingFrames.   
    """

    def __init__(self, inString="", minGene=0, largeOrf=False, startCodons={'ATG'}, stopCodons={'TAA', 'TAG', 'TGA'}):
        """
        Creates and initializes an OrfFinder by finding ORFs and adding them to openReadingFrames
        
        openReadingFrames contained all the information needed for the formatted output file so
        that was the only class attribute made and is accessible throughout the class.
        
        The constructor uses the methods, reverseComplimentaryStrand() and findOrfs() to add ORFs
        to the openReadingFrames list that will be used to generate a formatted output file.
        
        Keyword Arguments:
        self -- OrfFinder object
        inString -- A string of a DNA sequence (default - empty string)
        minGene -- Minimum length of a gene (default - 100)
        largeOrf -- Decides if the largest ORF or all ORFs are added (default -- False, Only largest ORF reported)
        startCodons -- Start codons used to find starts of ORFs (default - {'ATG'})
        stopCodons -- Stop codons used to find stops of ORFS (default - {'TAA', 'TAG', 'TGA'})
        """

        # Creates a complimentary strand of DNA
        complimentStrand = self.reverseComplimentaryStrand(inString)
        self.openReadingFrames = []
        # Finds and adds ORFs for DNA sequence
        self.findOrfs(inString, 1, minGene, startCodons, stopCodons, largeOrf)
        # Finds and adds ORFs for reverse DNA sequence
        self.findOrfs(complimentStrand, -1, minGene, startCodons, stopCodons, largeOrf)

    def reverse_comp_strand(self, fwd_strand):
        """ Creates a reversed complimentary strand of a DNA sequence. """
        complimentStrand = fwd_strand.lower().replace('a', 'T').replace('t', 'A').replace('c', 'G').replace('g',
                                                                                                            'C')
        return complimentStrand[::-1]  # Reverses and returns the new complimentary strand of DNA

    def findOrfs(self, strand, typeStrand, minLength, startCodons, stopCodons, largeOrf):
        """
        Generates ORFs that are added to openReadingFrames
        
        findOrfs() mainly decides what type of ORFs will be added to the list of ORFS.
        It uses a few for loops amd many parameters to determine start and stop position
        of an ORF based on the user's specifications.
        
        input: A string of a DNA sequence that will be searched for ORFs
        output: None
        
        The method creates a list of values that represent an ORF. It consists of
        Frame -- the frame the ORF was found
        Start Position -- The position the ORF begins on (including the Start Codon(s))
        Stop Position -- The position the ORF ends on (including the Stop Codon(s))
        Length -- The length of the ORF
        
        Keyword Arguments:
        Most of the keyword arguments of these methods were passed in the constructor.
        These arguments are the same as the same ones passed in.
        typeStrand -- an integer that determines if it a DNA sequence is forward DNA
        or the reversed DNA sequence.
        """
        # Loops through frames of a DNA sequence
        for frame in range(1, 4, 1):
            startPos = []  # List of start positions on the DNA sequence
            noStopYet = True  # Determines if stop codon has been reached
            # Used for boundary stop codon cases that occur in the beginning
            # Loops through codons
            for position in range(frame - 1, len(strand), 3):
                # Checks for start codons and if they can be added based on largeOrf argument
                if strand[position:position + 3] in startCodons and (largeOrf == False or len(startPos) == 0):
                    startPos.append(position + 1)  # Adds start positions to startPos[]
                # Checks for stop codons
                elif strand[position:position + 3] in stopCodons:
                    endPos = position + 3  # Finds end position
                    # Checks for boundary stop case ORFs
                    if noStopYet == True and len(startPos) == 0 and minLength < endPos:
                        if typeStrand > 0:  # Checks for positive frames
                            self.openReadingFrames.append([frame * typeStrand, 1, endPos, endPos])
                        elif typeStrand < 0:  # Checks for negative frames
                            self.openReadingFrames.append(
                                [frame * typeStrand, len(strand) - (endPos - 1), len(strand), endPos])
                    # Checks if start codons/positions upstream after previous stop
                    elif len(startPos) != 0:
                        for pos in startPos:
                            # Ensures ORFs are larger than minLength
                            if minLength < (endPos - pos + 1):
                                if typeStrand > 0:  # Checks for positive frames
                                    self.openReadingFrames.append([frame * typeStrand, pos, endPos, endPos - pos + 1])
                                else:
                                    self.openReadingFrames.append(
                                        [frame * typeStrand, len(strand) - (endPos - 1), len(strand) - (pos - 1),
                                         endPos - pos + 1])
                    noStopYet = False  # Does not allow for more boundary stop cases
                    startPos = []  # Gets rid of start positions after a stop is reached
            # Checks for boundary start cases
            if len(startPos) != 0:
                for pos in startPos:
                    # Ensures boundary start ORF cases are larger than minLength
                    if minLength < (len(strand) - frame + 2 - pos):
                        if typeStrand > 0:  # Checks for positive frames
                            self.openReadingFrames.append([frame * typeStrand, pos, len(strand), len(strand) - pos + 1])
                        else:
                            self.openReadingFrames.append(
                                [frame * typeStrand, 1, len(strand) - (pos - 1), len(strand) - pos + 1])
            # Check if no genes were found
            elif noStopYet == True and minLength < len(strand):
                # Adds whole DNA sequence as gene if no genes found
                self.openReadingFrames.append([frame * typeStrand, 1, len(strand), len(strand)])

    def getOrfs(self):
        """
        Returns ordered list of ORFs found on the DNA sequence
        
        input: None
        output: A list of ORFs that contain their start and stop positions, their frame, and their length
        
        The method sorts the ORFs by
        1. length: descending order
        2. frame: ascending order
        3. start position: ascending order
        """
        self.openReadingFrames.sort(key=lambda entry: (entry[3], -entry[0], -entry[1]), reverse=True)
        return self.openReadingFrames
