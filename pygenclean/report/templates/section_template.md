{{ "#" * section_level }} {{ section_name }}{{ " (" + reason + ")" if reason }} {{ "{#" }}sec-{{ section_key }}-{{ label_prefix }}}

{% if parent_step == "0" %}
This step uses the initial binary _Plink_ files.
{% else %}
The binary _Plink_ files used come from step {{ parent_step }} (see
@{{ parent_link }} for more information).
{% endif %}

{% block section_content %}{% endblock %}
