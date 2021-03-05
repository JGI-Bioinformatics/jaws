import os
import configparser


class Configuration():
    def __init__(self, config_file):
        """Read the given config file and create a config parser object containing the
        contents of the config file.

        :param config_file: file path and name of config file
        :type config_file: str
        :return: none
        """

        if not config_file:
            config_file = os.environ.get('JAWS_MONITOR_CONFIG')
        if not os.path.isfile(config_file):
            raise

        self.config = configparser.ConfigParser()
        try:
            self.config.read(config_file)
        except Exception:
            raise

    def get_port(self):
        """Parse the config entries to get the port number for serving up this web server using
        the prometheus client module.

        :param none
        :type none
        :return: none
        """

        port = self.config['MONITOR'].get('port')
        if port:
            return int(port)
        else:
            return None

    def get_config_section(self, section_name, is_partial_name=False):
        """Parse a config section and yield each entry under that section if is_partial_name=False.
        If is_partial_name=True, look for section that starts with the input section name and return
        a dictionary containing the entries under that section.

        :param section_name: name of section in config file
        :type section_name: string
        :param is_partial_name: boolean to check if input section name is a partial name. If False,
            yield all entries one at a time for each entry under the section name. If True,
            yield all entries one at a time in which the config section starts with the input
            section name. The entry returned is a dictionary containing the key/value pairs of
            the entries under that section.
        :type is_partial_name: boolean
        :return name: if is_partial_name=False, name of entry under the config section. If
            is_partial_name=True, name of section in lowercase stripping away the input section_name.
        :rtype name: string
        :return value: value of entry in config file.
        :rtype value: if is_partial_name=False, string. If_partial_name=True, dictionary.
        """

        if is_partial_name:
            for section in self.config.sections():
                if not section.startswith(section_name):
                    continue
                metric_name = section.replace(section_name, '', 1).lower()
                entries = dict(self.config.items(section))
                yield metric_name, entries
        else:
            for name, value in self.config.items(section_name):
                yield name, value
