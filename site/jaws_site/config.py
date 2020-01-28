"""
JAWS configuration class.  It imposes some restrictions upon the INI files, but facilitates creation, validation, and updates.  Namely, all sections and options must be defined in the template; no variable section or option names are allowed.
"""

import configparser
import sys
import os

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

    def validate(self, template_file):
        """
        Validate config file by comparing to a template file.  Sections and keys are added or removed as needed.
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
        changed = False

        # FIND AND ADD ANY MISSING SECTIONS
        section_printed = None
        for section in template.sections():
            section_printed = False
            if section not in self.sections():
                print("** ADDING TO CONFIG:\n[%s]" % (section,))
                changed = True
                section_printed = True
                self[section] = {}
            for key in template[section]:
                if key not in self[section]:
                    if not section_printed:
                        print("** ADDING TO CONFIG:\n[%s]" % (section,))
                        section_printed = True
                    print("%s = %s" % (key, template[section][key]))
                    changed = True
                    self[section][key] = template[section][key]
            
        # FIND AND DELETE ANY EXTRA/DEPRECATED SECTIONS
        for section in self.sections():
            if section not in template.sections():
                print("** DELETING FROM CONFIG:\n[%s]" % (section,))
                for key in self[section]:
                    print("%s = %s" % (key, self[section][key]))
                changed = True
                self.remove_section(section)
        for section in self.sections():
            section_printed = False
            for key in self[section]:
                if key not in template[section]:
                    if not section_printed:
                        print("** DELETING FROM CONFIG:\n[%s]" % (section, ))
                        section_printed = True
                    print("%s = %s" % (key, self[section][key]))
                    changed = True
                    self.remove_option(section, key)

        if changed:
            print("Writing new config file")
            with open(self.config_file, 'w') as f:
                self.write(f)
            sys.exit("Please review changes")

