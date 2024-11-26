import sys


class OpticConfigFile:
	"""
	A class to handle reading and validating a OPTIC configuration file.
	"""
	
	def __init__(self, config_file):
		"""
		Initializes the ConfigReader with the given configuration file.

		Args:
			config_file (str): Path to the configuration file.
		"""
		self.config_file = config_file
		self.attributes = {}
		self.config_options = {
			'File Extension',
			'Gene',
			'Start_position',
			'Reference Allele',
			'Alternate Allele',
			'Variant CDS',
			'Variant Amino Acid',
			'Variant Type',
			'Filter Column 1',
			'Filter Column 1 Exclusions',
			'Filter Column 2',
			'Filter Column 2 Exclusions'
		}
		self.int_keys = {'Gene',
		                 'Start_position',
		                 'Reference Allele',
		                 'Alternate Allele',
		                 'Variant CDS',
		                 'Variant Amino Acid',
		                 'Variant Type',
		                 'Filter Column 1',
		                 'Filter Column 2'}
		
		self.str_keys = {'File Extension',
		                 'Filter Column 1 Exclusions',
		                 'Filter Column 2 Exclusions'}
	
	def read_config(self):
		"""
		Reads and parses the configuration file, storing attributes.
		"""
		with open(self.config_file) as f:
			for line in f:
				if line.startswith('#') or not line.rstrip().split():
					continue
				line = line.rstrip().split('\t')
				self.attributes[line[0]] = line[1]
	
	def validate_config(self, args):
		"""
		Validates the configuration file for missing attributes and ensures correct types.

		Raises:
			SystemExit: If any required attributes are missing.
		"""
		
		in_config = set(self.attributes.keys())
		missing_keys = self.config_options - in_config
		errors = []
		
		# Check for specific missing keys
		for missing in missing_keys:
			if missing in ['File Extension', 'Gene', 'Start_position', 'Reference Allele', 'Alternate Allele']:
				errors.append(f"Missing required attribute in config: {missing}")
		
		# Check for one of either `Variant CDS` or `Variant Amino Acid`
		if 'Variant CDS' in missing_keys and 'Variant Amino Acid' in missing_keys:
			errors.append("One of either `Variant CDS` or `Variant Amino Acid` is required")
		
		if args.use_cds and 'Variant CDS' in missing_keys:
			errors.append("Cannot use `--use_cds` if no CDS value is provided in the config file")
		
		# Check dependencies for 'Filter Column 1 Exclusions'
		if 'Filter Column 1 Exclusions' in missing_keys and 'Filter Column 1' not in missing_keys:
			errors.append("'Filter Column 1 Exclusions' cannot be used without 'Filter Column 1' being set")
		
		# Check dependencies for 'Filter Column 2 Exclusions'
		if 'Filter Column 2 Exclusions' in missing_keys and 'Filter Column 2' not in missing_keys:
			errors.append("'Filter Column 2 Exclusions' cannot be used without 'Filter Column 2' being set")
		
		# Convert numeric strings to integers
		for key in self.int_keys:
			if key in self.attributes:
				try:
					self.attributes[key] = int(self.attributes[key])  # Convert to int
				except ValueError:
					errors.append(f"Key '{key}' must have a value of type int. Found: {self.attributes[key]}")
		
		# Convert string columns to string
		for key in self.str_keys:
			if key in self.attributes:
				self.attributes[key] = str(self.attributes[key])  # Convert to str
		
		if errors:
			for error in errors:
				print(error)
			sys.exit(1)
	
	def get_attributes(self):
		"""
		Returns the processed configuration attributes.

		Returns:
			dict: The configuration attributes.
		"""
		return self.attributes
