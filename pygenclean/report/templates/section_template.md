{{ "#" * section_level }} {{ section_name }}{{ " (" + reason + ")" if reason }} {{ "{#" }}sec-{{ section_key }}-{{ label_prefix }}}

{% block section_content %}{% endblock %}
