from src_edtop.models.atomic_model_cp import AtomicModelCP


def test_guiopengl(guiopengl_widget):
    widget = guiopengl_widget
    assert type(widget.main_model) == AtomicModelCP
