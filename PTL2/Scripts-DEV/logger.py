

class Logger():

	OKBLUE = '\033[94m'
	OKCYAN = '\033[96m'
	OKGREEN = '\033[92m'
	YELLOW = '\033[93m'
	RED = '\033[91m'
	ENDC = '\033[0m'
	BOLD = '\033[1m'
	UNDERLINE = '\033[4m'

	@staticmethod
	def Error(text, style = ''):

		print(Logger.RED + Logger.BOLD + "ERROR: " + Logger.ENDC + Logger.RED + style + text + Logger.ENDC + '\n')
		quit()

	@staticmethod
	def Warning(text, style = ''):

		print(Logger.YELLOW + Logger.BOLD + "WARNING: " + Logger.ENDC + Logger.YELLOW + style + text + Logger.ENDC + '\n')

	@staticmethod
	def Message(text, style = ''):

		print(Logger.OKGREEN + style + text + Logger.ENDC + '\n')
