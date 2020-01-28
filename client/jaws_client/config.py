"""
JAWS configuration class.  It imposes some restrictions upon the INI files, but facilitates creation, validation, and updates.  Namely, all sections and options must be defined in the template; no variable section or option names are allowed.
"""

import configparser
import sys
import os
import stat

DEBUG = False
if "JAWS_DEBUG" in os.environ: DEBUG = True

class Config(configparser.ConfigParser):

    def __init__(self, path=None, env=None, template=None):
        """
        Read config file or create new.  Validate parameters.  Update if necessary.
        """
        configparser.ConfigParser.__init__(self)
        if path:
            self.config_file = os.path.abspath(path)
        elif env:
            if env not in os.environ:
                sys.exit('Env var "%s" not defined' % (env,))
            self.config_file = os.environ[env]
        else:
            sys.exit('Either "path" or "env" required')

        if not os.path.isfile(self.config_file):
            self.create()
        try:
            self.read(self.config_file)
        except:
            sys.exit("Invalid config file: %s" % (self.config_file,))
        self.validate(template)
        

    def create(self):
        """
        Create a new config file
        """
        print("Creating new config file: %s" % (self.config_file,))
        with open(self.config_file, 'w') as f:
            self.write(f)
        os.chmod(self.config_file, stat.S_IRUSR|stat.S_IWUSR)


    def validate(self, template_file):
        """
        Validate config file by comparing to a template file.
        Sections and keys are added or removed as needed; this makes is easier to maintain changing config files in the wild.
        """
        # CHECK PACKAGE DIR FOR TEMPLATE "config.ini"
        if not template_file:
            template_file = "config.ini"
        if template_file == os.path.basename(template_file):
            template_file = os.path.join(os.path.abspath(os.path.dirname(__file__)), template_file)
        if not os.path.exists(template_file):
            sys.exit("Template file not found: %s" % (template_file,))
        template = configparser.ConfigParser()
        template.read(template_file)

        # FIND AND ADD ANY MISSING SECTIONS
        new_parameters = []
        for section in template.sections():
            if section not in self.sections():
                self[section] = {}
            for key in template[section]:
                if key not in self[section]:
                    new_parameters.append([section, key, template[section][key]])
                    self[section][key] = template[section][key]
        if new_parameters and DEBUG:
            print("** ATTENTION ** New configuration parameters were added (you may need to edit the file to set their values):")
            for section, key, value in new_parameters:
                print("[%s] %s = %s" % (section, key, value))
            
        # FIND AND DELETE ANY EXTRA/DEPRECATED SECTIONS
        removed_parameters = []
        for section in self.sections():
            if section not in template.sections():
                for key in self[section]:
                    removed_parameters.append([section, key, self[section][key]])
                self.remove_section(section)
        for section in self.sections():
            for key in self[section]:
                if key not in template[section]:
                    removed_parameters.append([section, key, self[section][key]])
                    self.remove_option(section, key)
        if removed_parameters and DEBUG:
            print("** ATTENTION ** Deprecated configuration parameters were removed:")
            for section, key, value in removed_parameters:
                print("[%s] %s = %s" % (section, key, value))

        if new_parameters or removed_parameters:
            print("Updating user config file at %s" % (self.config_file,))
            self.write_new_config_file()
            if DEBUG: sys.exit("Please review changes before retrying your command.")


    def write_new_config_file(self):
        """
        Write configuration parameters to file (i.e. after making changes).
        """
        with open(self.config_file, 'w') as f:
            self.write(f)
